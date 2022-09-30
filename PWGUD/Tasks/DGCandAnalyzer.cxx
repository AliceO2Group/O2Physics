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
// \brief Analyses reduced tables (DGCandidates, DGTracks) of DG candidates produced with DGCandProducer
//
//     options:
//           anaPars.mNCombine(2)
//           anaPars.mTPCnSigmas(120, 0.)
//
//           mTPCnSigmas contains 10 blocks (particles) of 12 elements:
//              0: PID
//              1: sign
//           2, 3: min/max nsigma for e
//           4, 5: min/max nsigma for pi
//           6, 7: min/max nsigma for mu
//           8, 9: min/max nsigma for Ka
//          10,11: min/max nsigma for Pr
//          In test for particle with PID it is required: min < nsigma < max
//          In test for all other particles it is required: nsigam < min || nsigam > max
//
//     usage: copts="--configuration json://DGCandAnalyzerConfig.json -b"
//
//           o2-analysis-ud-dgcand-analyzer $copts > DGCandAnalyzer.log
//
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  06.06.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Common/DataModel/PIDResponse.h"
#include "EventFiltering/PWGUD/DGCutparHolder.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/DGPIDSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // get a DGCutparHolder and DGAnaparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  DGAnaparHolder anaPars = DGAnaparHolder();
  MutableConfigurable<DGAnaparHolder> DGPars{"AnaPars", {}, "Analysis parameters"};

  ConfigurableAxis IVMAxis{"IVMAxis", {350, 0.0, 3.5}, ""};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, ""};

  // PID selector
  DGPIDSelector pidsel = DGPIDSelector();

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"nIVMs", "#nIVMs", {HistType::kTH1F, {{36, -0.5, 35.5}}}},
      {"TPCsignal1", "#TPCsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 100.0}}}},
      {"TPCsignal2", "#TPCsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 100.0}}}},
      {"sig1VsSig2TPC", "#sig1VsSig2TPC", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}}},
      {"TOFsignal1", "#TOFsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 50000.0}}}},
      {"TOFsignal2", "#TOFsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 50000.0}}}},
      {"sig1VsSig2TOF", "#sig1VsSig2TOF", {HistType::kTH2F, {{100, 0., 50000.}, {100, 0., 50000.}}}},
      {"nSigmaTPCPtEl", "#nSigmaTPCPtEl", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPi", "#nSigmaTPCPtPi", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtMu", "#nSigmaTPCPtMu", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtKa", "#nSigmaTPCPtKa", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPr", "#nSigmaTPCPtPr", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTOFPtEl", "#nSigmaTOFPtEl", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTOFPtPi", "#nSigmaTOFPtPi", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTOFPtMu", "#nSigmaTOFPtMu", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTOFPtKa", "#nSigmaTOFPtKa", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTOFPtPr", "#nSigmaTOFPtPr", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
    }};

  void fillNSigmaHist(DGParticle ivm, UDTracksFull dgtracks, DGPIDSelector pidsel)
  {
    // process only events with 2 tracks
    if (ivm.trkinds().size() != 2) {
      return;
    }

    // get analysis paraneters
    auto anaPars = pidsel.getAnaPars();

    // fill histogram
    auto tr1 = dgtracks.rawIteratorAt(ivm.trkinds()[0]);
    auto signalTPC1 = tr1.tpcSignal();
    auto tr2 = dgtracks.rawIteratorAt(ivm.trkinds()[1]);
    auto signalTPC2 = tr2.tpcSignal();

    registry.get<TH2>(HIST("TPCsignal1"))->Fill(tr1.pt(), signalTPC1);
    registry.get<TH2>(HIST("TPCsignal2"))->Fill(tr2.pt(), signalTPC2);
    registry.get<TH2>(HIST("sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);

    auto signalTOF1 = tr1.tofSignal();
    auto signalTOF2 = tr2.tofSignal();

    registry.get<TH2>(HIST("TOFsignal1"))->Fill(tr1.pt(), signalTOF1);
    registry.get<TH2>(HIST("TOFsignal2"))->Fill(tr2.pt(), signalTOF2);
    registry.get<TH2>(HIST("sig1VsSig2TOF"))->Fill(signalTOF1, signalTOF2);
  }

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
    anaPars = (DGAnaparHolder)DGPars;
    pidsel.init(anaPars);
    if (verbose) {
      pidsel.Print();
    }

    const AxisSpec axisIVM{IVMAxis, "IVM axis for histograms"};
    const AxisSpec axispt{ptAxis, "pt axis for histograms"};
    registry.add("IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {axisIVM, axispt}});
    registry.add("IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {axisIVM, axispt}});
  }

  void process(aod::UDCollision const& dgcand, UDTracksFull const& dgtracks)
  {

    // skip events with too few/many tracks
    if (dgcand.numContrib() < diffCuts.minNTracks() || dgcand.numContrib() > diffCuts.maxNTracks()) {
      return;
    }

    // skip events with out-of-range net charge
    if (dgcand.netCharge() < diffCuts.minNetCharge() || dgcand.netCharge() > diffCuts.maxNetCharge()) {
      return;
    }

    // skip events with out-of-range rgtrwTOF
    if (dgcand.rgtrwTOF() < diffCuts.minRgtrwTOF()) {
      return;
    }

    // find track combinations which are compatible with PID cuts
    auto nIVMs = pidsel.computeIVMs(dgtracks);

    // update histograms
    registry.get<TH1>(HIST("nIVMs"))->Fill(nIVMs, 1.);
    for (auto ivm : pidsel.IVMs()) {
      registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = dgtracks.rawIteratorAt(ind);
        registry.get<TH2>(HIST("IVMptTrkDG"))->Fill(ivm.M(), track.pt());

        // fill nSigma histograms
        registry.get<TH2>(HIST("nSigmaTPCPtEl"))->Fill(track.pt(), track.tpcNSigmaEl());
        registry.get<TH2>(HIST("nSigmaTPCPtPi"))->Fill(track.pt(), track.tpcNSigmaPi());
        registry.get<TH2>(HIST("nSigmaTPCPtMu"))->Fill(track.pt(), track.tpcNSigmaMu());
        registry.get<TH2>(HIST("nSigmaTPCPtKa"))->Fill(track.pt(), track.tpcNSigmaKa());
        registry.get<TH2>(HIST("nSigmaTPCPtPr"))->Fill(track.pt(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          registry.get<TH2>(HIST("nSigmaTOFPtEl"))->Fill(track.pt(), track.tofNSigmaEl());
          registry.get<TH2>(HIST("nSigmaTOFPtPi"))->Fill(track.pt(), track.tofNSigmaPi());
          registry.get<TH2>(HIST("nSigmaTOFPtMu"))->Fill(track.pt(), track.tofNSigmaMu());
          registry.get<TH2>(HIST("nSigmaTOFPtKa"))->Fill(track.pt(), track.tofNSigmaKa());
          registry.get<TH2>(HIST("nSigmaTOFPtPr"))->Fill(track.pt(), track.tofNSigmaPr());
        }
      }
      fillNSigmaHist(ivm, dgtracks, pidsel);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandAnalyzer>(cfgc, TaskName{"dgcandanalyzer"}),
  };
}
