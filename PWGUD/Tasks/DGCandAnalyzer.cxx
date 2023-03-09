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
#include "PWGUD/Core/UDFSParser.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  Configurable<int> candCaseSel{"CandCase", {}, "0: all Cands, 1: only ColCands,2: only BCCands"};
  Configurable<std::string> goodRunsFile{"goodRunsFile", {}, "json with list of good runs"};
  Configurable<std::string> fillingSchemeFile{"fillingSchemeFile", {}, "csv file with filling scheme information"};

  // a pdg object
  TDatabasePDG* pdg = nullptr;

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

  // a filling scheme parser
  UDFSParser fsparser = UDFSParser();

  // a global container to contain bcnum of accepted candidates
  std::set<uint64_t> bcnums;

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"nIVMs", "#nIVMs", {HistType::kTH1F, {{36, -0.5, 35.5}}}},
      {"candCase", "#candCase", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
      {"TPCsignal1", "#TPCsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, -100., 100.0}}}},
      {"TPCsignal2", "#TPCsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, -100., 100.0}}}},
      {"sig1VsSig2TPC", "#sig1VsSig2TPC", {HistType::kTH2F, {{100, -100., 100.}, {100, -100., 100.}}}},
      {"TOFsignal1", "#TOFsignal1", {HistType::kTH2F, {{100, 0., 3.}, {400, -200., 200.}}}},
      {"TOFsignal2", "#TOFsignal2", {HistType::kTH2F, {{100, 0., 3.}, {400, -200., 200.}}}},
      {"sig1VsSig2TOF", "#sig1VsSig2TOF", {HistType::kTH2F, {{100, -200., 200.}, {100, -200., 200.}}}},
      {"nSigmaTPCPtEl", "#nSigmaTPCPtEl", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPi", "#nSigmaTPCPtPi", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtMu", "#nSigmaTPCPtMu", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtKa", "#nSigmaTPCPtKa", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"nSigmaTPCPtPr", "#nSigmaTPCPtPr", {HistType::kTH2F, {{250, 0.0, 2.5}, {100, -20.0, 20.0}}}},
      {"eta1Vseta2", "#eta1Vseta2", {HistType::kTH2F, {{200, -2.0, 2.0}, {200, -2.0, 2.0}}}},
    }};

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  template <typename TTrack>
  void fillSignalHists(DGParticle ivm, TTrack const& dgtracks, DGPIDSelector pidsel)
  {
    // process only events with 2 tracks
    if (ivm.trkinds().size() != 2) {
      return;
    }

    // fill histogram
    auto tr1 = dgtracks.begin() + ivm.trkinds()[0];
    auto signalTPC1 = tr1.tpcSignal() * tr1.sign();
    auto tr2 = dgtracks.begin() + ivm.trkinds()[1];
    auto signalTPC2 = tr2.tpcSignal() * tr2.sign();

    registry.get<TH2>(HIST("TPCsignal1"))->Fill(tr1.pt(), signalTPC1);
    registry.get<TH2>(HIST("TPCsignal2"))->Fill(tr2.pt(), signalTPC2);
    registry.get<TH2>(HIST("sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);

    auto signalTOF1 = tr1.tofSignal() / 1.E3;
    auto signalTOF2 = tr2.tofSignal() / 1.E3;

    registry.get<TH2>(HIST("TOFsignal1"))->Fill(tr1.pt(), signalTOF1);
    registry.get<TH2>(HIST("TOFsignal2"))->Fill(tr2.pt(), signalTOF2);
    registry.get<TH2>(HIST("sig1VsSig2TOF"))->Fill(signalTOF1, signalTOF2);

    auto m1 = particleMass(pdg, pidsel.getAnaPars().PIDs()[0]);
    auto ene1 = sqrt(pow(tr1.px(), 2.) + pow(tr1.py(), 2.) + pow(tr1.pz(), 2.) + m1);
    auto lv1 = TLorentzVector(tr1.px(), tr1.py(), tr1.pz(), ene1);
    LOGF(debug, "pid1 %f mass %f energy %f", pidsel.getAnaPars().PIDs()[0], m1, ene1);
    auto m2 = particleMass(pdg, pidsel.getAnaPars().PIDs()[1]);
    auto ene2 = sqrt(pow(tr2.px(), 2.) + pow(tr2.py(), 2.) + pow(tr2.pz(), 2.) + m2);
    auto lv2 = TLorentzVector(tr2.px(), tr2.py(), tr2.pz(), ene2);
    LOGF(debug, "pid2 %f mass %f energy %f", pidsel.getAnaPars().PIDs()[1], m2, ene2);
    registry.get<TH2>(HIST("eta1Vseta2"))->Fill(lv1.Eta(), lv2.Eta());
  }

  void init(InitContext&)
  {
    pdg = TDatabasePDG::Instance();

    diffCuts = (DGCutparHolder)DGCuts;
    anaPars = (DGAnaparHolder)DGPars;
    pidsel.init(anaPars);
    grsel.init(goodRunsFile);
    std::string FSFile(fillingSchemeFile);
    fsparser.readFS(FSFile.data());

    if (verbose) {
      pidsel.Print();
      grsel.Print();
      fsparser.Print();
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
    registry.add("2Tracketa1IVM", "#2Tracketa1IVM", {HistType::kTH2F, {axisIVM, {160, -2.0, 2.0}}});
    registry.add("2Tracketa2IVM", "#2Tracketa2IVM", {HistType::kTH2F, {axisIVM, {160, -2.0, 2.0}}});

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

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    for (auto track : dgtracks) {
      auto mom = sqrt(pow(track.pt(), 2) + pow(track.pz(), 2));
      LOGF(debug, "<DGCandAnalyzer> %d %f %f %f %f %f", track.isPVContributor(), track.px(), track.py(), track.pz(), track.pt(), mom);
    }
    LOGF(debug, "");

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
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);
    if (dgcand.numContrib() != PVContributors.size()) {
      LOGF(info, "Missmatch of PVContributors %d != %d", dgcand.numContrib(), PVContributors.size());
    }
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
    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    auto minRgtrwTOF = candCase != 1 ? 1.0 : diffCuts.minRgtrwTOF();
    if (rtrwTOF < minRgtrwTOF) {
      LOGF(debug, "Rejected 3: %f below threshold of %f.", rtrwTOF, minRgtrwTOF);
      return;
    }

    // check FIT information
    auto bitMin = anaPars.dBCMin() + 16;
    auto bitMax = anaPars.dBCMax() + 16;
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (TESTBIT(dgcand.bbFT0Apf(), bit) ||
          TESTBIT(dgcand.bbFT0Cpf(), bit) ||
          TESTBIT(dgcand.bbFV0Apf(), bit) ||
          TESTBIT(dgcand.bbFDDApf(), bit) ||
          TESTBIT(dgcand.bbFDDCpf(), bit)) {
        return;
      }
    }

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
    auto nIVMs = pidsel.computeIVMs(PVContributors);

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

      // is BB bunch?
      if (!fsparser.isP2BCBB(bcnum % o2::constants::lhc::LHCMaxBunches)) {
        LOGF(info, "bcnum %d is not a BB BC", bcnum % o2::constants::lhc::LHCMaxBunches);
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
      if (anaPars.nCombine() == 2) {
        auto ind1 = ivm.trkinds()[0];
        auto trk1 = PVContributors.begin() + ind1;
        auto v1 = TVector3(trk1.px(), trk1.py(), trk1.pz());
        auto ind2 = ivm.trkinds()[1];
        auto trk2 = PVContributors.begin() + ind2;
        auto v2 = TVector3(trk2.px(), trk2.py(), trk2.pz());

        // cut on angle
        auto angle = v1.Angle(v2);
        LOGF(debug, "angle %f (%f / %f)", angle, anaPars.minAlpha(), anaPars.maxAlpha());
        if (angle < anaPars.minAlpha() || angle > anaPars.maxAlpha()) {
          continue;
        } else {
          registry.get<TH1>(HIST("2TrackAngle"))->Fill(angle, 1.);
          registry.get<TH2>(HIST("2TrackAngleIVM"))->Fill(ivm.M(), angle, 1.);
        }

        registry.get<TH2>(HIST("2Tracketa1IVM"))->Fill(ivm.M(), v1.Eta(), 1.);
        registry.get<TH2>(HIST("2Tracketa2IVM"))->Fill(ivm.M(), v2.Eta(), 1.);
      }

      registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = PVContributors.begin() + ind;
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
      // fillSignalHists(ivm, dgtracks, pidsel);
      fillSignalHists(ivm, PVContributors, pidsel);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandAnalyzer>(cfgc, TaskName{"dgcandanalyzer"}),
  };
}
