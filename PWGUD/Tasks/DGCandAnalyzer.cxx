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
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/DGCutparHolder.h"
#include "PWGUD/Core/DGPIDSelector.h"
#include "PWGUD/Core/UDGoodRunSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  Configurable<int> candCaseSel{"CandCase", {}, "0: all Cands, 1: only ColCands,2: only BCCands"};
  Configurable<std::string> goodRunsFile{"goodRunsFile", {}, "json with list of good runs"};

  // get a DGCutparHolder and DGAnaparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // analysis cuts
  DGAnaparHolder anaPars = DGAnaparHolder();
  Configurable<DGAnaparHolder> DGPars{"anaPars", {}, "Analysis parameters"};

  ConfigurableAxis IVMAxis{"IVMAxis", {350, 0.0, 3.5}, ""};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, ""};
  ConfigurableAxis nsTOFAxis{"nsTOFAxis", {100, -100.0, 100.0}, ""};

  // PID and goodRun selector
  DGPIDSelector pidsel = DGPIDSelector();
  UDGoodRunSelector grsel = UDGoodRunSelector();

  // a global container to contain bcnum of accepted candidates
  std::set<uint64_t> bcnums;

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"nIVMs", "#nIVMs", {HistType::kTH1F, {{36, -0.5, 35.5}}}},
      {"candCase", "#candCase", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
      {"TPCsignal1", "#TPCsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 100.0}}}},
      {"TPCsignal2", "#TPCsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, 0., 100.0}}}},
      {"sig1VsSig2TPC", "#sig1VsSig2TPC", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}}},
      {"TOFsignal1", "#TOFsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, -1000., 1000.}}}},
      {"TOFsignal2", "#TOFsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, -1000., 1000.}}}},
      {"sig1VsSig2TOF", "#sig1VsSig2TOF", {HistType::kTH2F, {{100, -1000., 1000.}, {100, -1000., 1000.}}}},
      {"nSigmaTPCPtEl", "#nSigmaTPCPtEl", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPi", "#nSigmaTPCPtPi", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtMu", "#nSigmaTPCPtMu", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtKa", "#nSigmaTPCPtKa", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPr", "#nSigmaTPCPtPr", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
    }};

  void fillSignalHists(DGParticle ivm, UDTracksFull dgtracks, DGPIDSelector pidsel)
  {
    // process only events with 2 tracks
    if (ivm.trkinds().size() != 2) {
      return;
    }

    // fill histogram
    auto tr1 = dgtracks.rawIteratorAt(ivm.trkinds()[0]);
    auto signalTPC1 = tr1.tpcSignal();
    auto tr2 = dgtracks.rawIteratorAt(ivm.trkinds()[1]);
    auto signalTPC2 = tr2.tpcSignal();

    registry.get<TH2>(HIST("TPCsignal1"))->Fill(tr1.pt(), signalTPC1);
    registry.get<TH2>(HIST("TPCsignal2"))->Fill(tr2.pt(), signalTPC2);
    registry.get<TH2>(HIST("sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);

    auto signalTOF1 = tr1.tofSignal() / 1.E3;
    auto signalTOF2 = tr2.tofSignal() / 1.E3;

    registry.get<TH2>(HIST("TOFsignal1"))->Fill(tr1.pt(), signalTOF1);
    registry.get<TH2>(HIST("TOFsignal2"))->Fill(tr2.pt(), signalTOF2);
    registry.get<TH2>(HIST("sig1VsSig2TOF"))->Fill(signalTOF1, signalTOF2);
  }

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
    anaPars = (DGAnaparHolder)DGPars;
    pidsel.init(anaPars);
    grsel.init(goodRunsFile);

    if (verbose) {
      pidsel.Print();
      grsel.Print();
    }
    bcnums.clear();

    const AxisSpec axisIVM{IVMAxis, "IVM axis for histograms"};
    const AxisSpec axispt{ptAxis, "pt axis for histograms"};
    registry.add("trackQC", "#trackQC", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("dcaXYDG", "#dcaXYDG", {HistType::kTH1F, {{400, -2., 2.}}});
    registry.add("ptTrkdcaXYDG", "#ptTrkdcaXYDG", {HistType::kTH2F, {axispt, {80, -2., 2.}}});
    registry.add("dcaZDG", "#dcaZDG", {HistType::kTH1F, {{800, -20., 20.}}});
    registry.add("ptTrkdcaZDG", "#ptTrkdcaZDG", {HistType::kTH2F, {axispt, {400, -20., 20.}}});
    registry.add("IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {axisIVM, axispt}});
    registry.add("IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {axisIVM, axispt}});

    const AxisSpec axisnsTOF{nsTOFAxis, "nSigma TOF axis for histograms"};
    registry.add("nSigmaTOFPtEl", "#nSigmaTOFPtEl", {HistType::kTH2F, {{250, 0.0, 2.5}, axisnsTOF}});
    registry.add("nSigmaTOFPtPi", "#nSigmaTOFPtPi", {HistType::kTH2F, {{250, 0.0, 2.5}, axisnsTOF}});
    registry.add("nSigmaTOFPtMu", "#nSigmaTOFPtMu", {HistType::kTH2F, {{250, 0.0, 2.5}, axisnsTOF}});
    registry.add("nSigmaTOFPtKa", "#nSigmaTOFPtKa", {HistType::kTH2F, {{250, 0.0, 2.5}, axisnsTOF}});
    registry.add("nSigmaTOFPtPr", "#nSigmaTOFPtPr", {HistType::kTH2F, {{250, 0.0, 2.5}, axisnsTOF}});

    registry.add("2TrackAngle", "#2TrackAngle", {HistType::kTH1F, {{140, -0.2, 3.3}}});
    registry.add("2TrackAngleIVM", "#2TrackAngleIVM", {HistType::kTH2F, {axisIVM, {140, -0.2, 3.3}}});

    registry.add("FT0AAmplitude", "#FT0AAmplitude", {HistType::kTH1F, {{5000, 0., 5000.}}});
    registry.add("FT0CAmplitude", "#FT0CAmplitude", {HistType::kTH1F, {{5000, 0., 5000.}}});
    registry.add("FV0AAmplitude", "#FV0AAmplitude", {HistType::kTH1F, {{5000, 0., 5000.}}});
    registry.add("FDDAAmplitude", "#FDDAAmplitude", {HistType::kTH1F, {{5000, 0., 5000.}}});
    registry.add("FDDCAmplitude", "#FDDCAmplitude", {HistType::kTH1F, {{5000, 0., 5000.}}});

    registry.add("BBT0A", "#BBT0A", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("BBT0C", "#BBT0C", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("BBV0A", "#BBV0A", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("BBFDDA", "#BBFDDA", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("BBFDDC", "#BBFDDC", {HistType::kTH1F, {{32, -16.5, 15.5}}});
  }

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    // accept only selected run numbers
    if (!grsel.isGoodRun(dgcand.runNumber())) {
      return;
    }
    LOGF(debug, "Run number %d", dgcand.runNumber());

    // skip unwanted cases
    // 0. all candidates
    // 1. candidate has associated BC and associated collision
    // 2. candidate has associated BC but no associated collision
    // 3. candidate has no associated BC
    int candCase = 1;
    if (dgcand.posX() == -1. && dgcand.posY() == 1. && dgcand.posZ() == -1.) {
      candCase = 2;
    } else if (dgcand.posX() == -2. && dgcand.posY() == 2. && dgcand.posZ() == -2.) {
      candCase = 3;
    }
    if (candCaseSel > 0 && candCase != candCaseSel) {
      return;
    }

    // skip events with too few/many tracks
    if (dgcand.numContrib() < diffCuts.minNTracks() || dgcand.numContrib() > diffCuts.maxNTracks()) {
      LOGF(debug, "Rejected 1: %d not in range [%d, %d].", dgcand.numContrib(), diffCuts.minNTracks(), diffCuts.maxNTracks());
      return;
    }

    // skip events with out-of-range net charge
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), dgcand.netCharge()) == netChargeValues.end()) {
      LOGF(debug, "Rejected 2: %d not in set.", dgcand.netCharge());
      return;
    }

    // skip events with out-of-range rgtrwTOF
    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, dgtracks.size());
    auto minRgtrwTOF = candCase != 1 ? 1.0 : diffCuts.minRgtrwTOF();
    if (rtrwTOF < minRgtrwTOF) {
      LOGF(debug, "Rejected 3: %f below threshold of %f.", rtrwTOF, minRgtrwTOF);
      return;
    }

    // check FIT information
    // for (auto bit = 15; bit <= 17; bit++) {
    //  if (TESTBIT(dgcand.bbFT0Apf(), bit) ||
    //      TESTBIT(dgcand.bbFT0Cpf(), bit) ||
    //      TESTBIT(dgcand.bbFV0Apf(), bit) ||
    //      TESTBIT(dgcand.bbFDDApf(), bit) ||
    //      TESTBIT(dgcand.bbFDDCpf(), bit)) {
    //    return;
    //  }
    //}

    // fill FIT amplitude histograms
    registry.get<TH1>(HIST("FT0AAmplitude"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
    registry.get<TH1>(HIST("FT0CAmplitude"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
    registry.get<TH1>(HIST("FV0AAmplitude"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
    registry.get<TH1>(HIST("FDDAAmplitude"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
    registry.get<TH1>(HIST("FDDCAmplitude"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);

    // fill BBFlag histograms
    for (auto bit = 0; bit < 33; bit++) {
      registry.get<TH1>(HIST("BBT0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Apf(), bit));
      registry.get<TH1>(HIST("BBT0C"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Cpf(), bit));
      registry.get<TH1>(HIST("BBV0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFV0Apf(), bit));
      registry.get<TH1>(HIST("BBFDDA"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDApf(), bit));
      registry.get<TH1>(HIST("BBFDDC"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDCpf(), bit));
    }

    // find track combinations which are compatible with PID cuts
    auto nIVMs = pidsel.computeIVMs(dgtracks);

    // update candCase histogram
    if (nIVMs > 0) {
      registry.get<TH1>(HIST("candCase"))->Fill(candCase, 1.);
      // check bcnum
      auto bcnum = dgcand.globalBC();
      if (bcnums.find(bcnum) != bcnums.end()) {
        LOGF(info, "candCase %i bcnum %i allready found! ", candCase, bcnum);
        registry.get<TH1>(HIST("candCase"))->Fill(4, 1.);
        return;
      } else {
        bcnums.insert(bcnum);
      }
    } else {
      LOGF(debug, "Rejected 4: no IVMs.");
    }

    // update histograms
    registry.get<TH1>(HIST("nIVMs"))->Fill(nIVMs, 1.);
    for (auto ivm : pidsel.IVMs()) {

      // cut on pt-system
      if (ivm.Perp() < anaPars.minptsys() || ivm.Perp() > anaPars.maxptsys()) {
        continue;
      }

      // applicable to 2-track events - cut on angle between two tracks
      if (dgcand.numContrib() == 2) {
        auto ind1 = ivm.trkinds()[0];
        auto trk1 = dgtracks.rawIteratorAt(ind1);
        auto v1 = TVector3(trk1.px(), trk1.py(), trk1.pz());
        auto ind2 = ivm.trkinds()[1];
        auto trk2 = dgtracks.rawIteratorAt(ind2);
        auto v2 = TVector3(trk2.px(), trk2.py(), trk2.pz());

        auto angle = v1.Angle(v2);
        LOGF(debug, "angle %f", angle);

        // cut on angle
        if (angle < anaPars.minAlpha() || angle > anaPars.maxAlpha()) {
          continue;
        } else {
          registry.get<TH1>(HIST("2TrackAngle"))->Fill(angle, 1.);
          registry.get<TH2>(HIST("2TrackAngleIVM"))->Fill(ivm.M(), angle, 1.);
        }
      }

      registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = dgtracks.rawIteratorAt(ind);
        registry.get<TH1>(HIST("trackQC"))->Fill(0., 1.);
        registry.get<TH1>(HIST("trackQC"))->Fill(1., track.hasITS() * 1.);
        registry.get<TH1>(HIST("trackQC"))->Fill(2., track.hasTPC() * 1.);
        registry.get<TH1>(HIST("trackQC"))->Fill(3., track.hasTRD() * 1.);
        registry.get<TH1>(HIST("trackQC"))->Fill(4., track.hasTOF() * 1.);
        // registry.get<TH1>(HIST("dcaXYDG"))->Fill(track.dcaXY());
        // registry.get<TH2>(HIST("ptTrkdcaXYDG"))->Fill(track.pt(), track.dcaXY());
        // registry.get<TH1>(HIST("dcaZDG"))->Fill(track.dcaZ());
        // registry.get<TH2>(HIST("ptTrkdcaZDG"))->Fill(track.pt(), track.dcaZ());

        registry.get<TH2>(HIST("IVMptTrkDG"))->Fill(ivm.M(), track.pt());

        // fill nSigma histograms
        registry.get<TH2>(HIST("nSigmaTPCPtEl"))->Fill(track.pt(), track.tpcNSigmaEl());
        registry.get<TH2>(HIST("nSigmaTPCPtPi"))->Fill(track.pt(), track.tpcNSigmaPi());
        registry.get<TH2>(HIST("nSigmaTPCPtMu"))->Fill(track.pt(), track.tpcNSigmaMu());
        registry.get<TH2>(HIST("nSigmaTPCPtKa"))->Fill(track.pt(), track.tpcNSigmaKa());
        registry.get<TH2>(HIST("nSigmaTPCPtPr"))->Fill(track.pt(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          LOGF(debug, "tofNSigmaPi %f", track.tofNSigmaPi());
          registry.get<TH2>(HIST("nSigmaTOFPtEl"))->Fill(track.pt(), track.tofNSigmaEl());
          registry.get<TH2>(HIST("nSigmaTOFPtPi"))->Fill(track.pt(), track.tofNSigmaPi());
          registry.get<TH2>(HIST("nSigmaTOFPtMu"))->Fill(track.pt(), track.tofNSigmaMu());
          registry.get<TH2>(HIST("nSigmaTOFPtKa"))->Fill(track.pt(), track.tofNSigmaKa());
          registry.get<TH2>(HIST("nSigmaTOFPtPr"))->Fill(track.pt(), track.tofNSigmaPr());
        }
      }
      fillSignalHists(ivm, dgtracks, pidsel);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandAnalyzer>(cfgc, TaskName{"dgcandanalyzer"}),
  };
}
