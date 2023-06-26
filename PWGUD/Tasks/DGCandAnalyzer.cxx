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

#include <set>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "CommonConstants/LHCConstants.h"
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

  // ccdb
  Service<o2::ccdb::BasicCCDBManager> ccdb;

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

  // filling scheme
  int lastRun = -1;                                          // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB; // bc pattern of colliding bunches

  // a global container to contain bcnum of accepted candidates
  std::set<uint64_t> bcnums;

  // define histograms
  HistogramRegistry registry{
    "registry",
    {{"nDGperRun", "Number of DG collisions per run", {HistType::kTH1D, {{1, 0, 1}}}},
     {"nIVMs", "Number of IVMs per DG collision", {HistType::kTH1F, {{36, -0.5, 35.5}}}},
     {"candCase", "#candCase", {HistType::kTH1F, {{5, -0.5, 4.5}}}}}};

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
    auto m1 = particleMass(pdg, pidsel.getAnaPars().PIDs()[0]);
    auto ene1 = sqrt(pow(tr1.px(), 2.) + pow(tr1.py(), 2.) + pow(tr1.pz(), 2.) + m1);
    auto lv1 = TLorentzVector(tr1.px(), tr1.py(), tr1.pz(), ene1);
    LOGF(debug, "pid1 %f mass %f energy %f", pidsel.getAnaPars().PIDs()[0], m1, ene1);
    auto signalTPC1 = tr1.tpcSignal();
    auto signalTOF1 = tr1.tofSignal() / 1.E6;

    auto tr2 = dgtracks.begin() + ivm.trkinds()[1];
    auto m2 = particleMass(pdg, pidsel.getAnaPars().PIDs()[1]);
    auto ene2 = sqrt(pow(tr2.px(), 2.) + pow(tr2.py(), 2.) + pow(tr2.pz(), 2.) + m2);
    auto lv2 = TLorentzVector(tr2.px(), tr2.py(), tr2.pz(), ene2);
    LOGF(debug, "pid2 %f mass %f energy %f", pidsel.getAnaPars().PIDs()[1], m2, ene2);
    auto signalTPC2 = tr2.tpcSignal();
    auto signalTOF2 = tr2.tofSignal() / 1.E6;

    registry.get<TH2>(HIST("TPCsignal1"))->Fill(tr1.tpcInnerParam(), signalTPC1);
    registry.get<TH2>(HIST("TPCsignal2"))->Fill(tr2.tpcInnerParam(), signalTPC2);
    registry.get<TH2>(HIST("sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);
    registry.get<TH2>(HIST("eta1Vseta2"))->Fill(lv1.Eta(), lv2.Eta());

    if (tr1.hasTOF()) {
      registry.get<TH2>(HIST("TOFsignal1"))->Fill(lv1.P(), signalTOF1);
    }
    if (tr2.hasTOF()) {
      registry.get<TH2>(HIST("TOFsignal2"))->Fill(lv2.P(), signalTOF2);
    }
    if (tr1.hasTOF() && tr2.hasTOF()) {
      registry.get<TH2>(HIST("sig1VsSig2TOF"))->Fill(signalTOF1, signalTOF2);
    }
  }

  void init(InitContext&)
  {
    // initalise ccdb
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // PDG
    pdg = TDatabasePDG::Instance();

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
    registry.add("nSigmaTPCPEl", "#nSigmaTPCPEl", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
    registry.add("nSigmaTPCPPi", "#nSigmaTPCPPi", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
    registry.add("nSigmaTPCPMu", "#nSigmaTPCPMu", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
    registry.add("nSigmaTPCPKa", "#nSigmaTPCPKa", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
    registry.add("nSigmaTPCPPr", "#nSigmaTPCPPr", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});

    registry.add("trackQC", "#trackQC", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("dcaXYDG", "#dcaXYDG", {HistType::kTH1F, {{400, -2., 2.}}});
    registry.add("ptTrkdcaXYDG", "#ptTrkdcaXYDG", {HistType::kTH2F, {axispt, {80, -2., 2.}}});
    registry.add("dcaZDG", "#dcaZDG", {HistType::kTH1F, {{800, -20., 20.}}});
    registry.add("ptTrkdcaZDG", "#ptTrkdcaZDG", {HistType::kTH2F, {axispt, {400, -20., 20.}}});
    registry.add("IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {axisIVM, axispt}});
    registry.add("IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {axisIVM, axispt}});

    const AxisSpec axisnsTOF{nsTOFAxis, "nSigma TOF axis for histograms"};
    registry.add("nSigmaTOFPEl", "#nSigmaTOFPEl", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("nSigmaTOFPPi", "#nSigmaTOFPPi", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("nSigmaTOFPMu", "#nSigmaTOFPMu", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("nSigmaTOFPKa", "#nSigmaTOFPKa", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("nSigmaTOFPPr", "#nSigmaTOFPPr", {HistType::kTH2F, {axispt, axisnsTOF}});

    // FIT signals
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

    // 2 track events
    registry.add("TPCChi2NCl1", "#TPCChi2NCl1", {HistType::kTH1F, {{200, 0., 50.}}});
    registry.add("TPCChi2NCl2", "#TPCChi2NCl2", {HistType::kTH1F, {{200, 0., 50.}}});
    registry.add("TPCsignal1", "#TPCsignal1", {HistType::kTH2F, {{1000, 0., 10.}, {500, 0., 500.}}});
    registry.add("TPCsignal2", "#TPCsignal2", {HistType::kTH2F, {{1000, 0., 10.}, {500, 0., 500.}}});
    registry.add("sig1VsSig2TPC", "#sig1VsSig2TPC", {HistType::kTH2F, {{300, 0., 300.}, {300, 0., 300.}}});
    registry.add("TOFsignal1", "#TOFsignal1", {HistType::kTH2F, {{1000, 0., 10.}, {400, -100., 100.}}});
    registry.add("TOFsignal2", "#TOFsignal2", {HistType::kTH2F, {{1000, 0., 10.}, {400, -100., 100.}}});
    registry.add("sig1VsSig2TOF", "#sig1VsSig2TOF", {HistType::kTH2F, {{160, -20., 60.}, {160, -20., 60.}}});
    registry.add("eta1Vseta2", "#eta1Vseta2", {HistType::kTH2F, {{200, -2.0, 2.0}, {200, -2.0, 2.0}}});
    registry.add("2Trackpt1pt2", "#2Trackpt1pt2", {HistType::kTH2F, {axispt, axispt}});
    registry.add("2Trackpt1eta1", "#2Trackpt1eta1", {HistType::kTH2F, {axispt, {200, -2.0, 2.0}}});
    registry.add("2Trackpt2eta2", "#2Trackpt2eta2", {HistType::kTH2F, {axispt, {200, -2.0, 2.0}}});
    registry.add("2TrackAngle", "#2TrackAngle", {HistType::kTH1F, {{140, -0.2, 3.3}}});
    registry.add("2TrackAngleIVM", "#2TrackAngleIVM", {HistType::kTH2F, {axisIVM, {140, -0.2, 3.3}}});
    registry.add("2Tracketa1IVM", "#2Tracketa1IVM", {HistType::kTH2F, {axisIVM, {160, -2.0, 2.0}}});
    registry.add("2Tracketa2IVM", "#2Tracketa2IVM", {HistType::kTH2F, {axisIVM, {160, -2.0, 2.0}}});
  }

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    // count collisions
    registry.get<TH1>(HIST("candCase"))->Fill(0., 1.);

    // accept only selected run numbers
    int run = dgcand.runNumber();
    if (!grsel.isGoodRun(run)) {
      return;
    }

    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= 500000) {
      LOGF(info, "Updating bcPattern %d ...", run);
      auto tss = ccdb->getRunDuration(run);
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tss.first);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();
      lastRun = run;
      LOGF(info, "done!");
    }

    // is BB bunch?
    auto bcnum = dgcand.globalBC();
    if (run >= 500000 && bcPatternB[bcnum % o2::constants::lhc::LHCMaxBunches] == 0) {
      LOGF(debug, "bcnum[1] %d is not a BB BC", bcnum % o2::constants::lhc::LHCMaxBunches);
      return;
    }

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
    LOGF(debug, "Number of IVMs %d", nIVMs);

    // update candCase histogram
    if (nIVMs > 0) {
      registry.get<TH1>(HIST("candCase"))->Fill(candCase, 1.);
      // check bcnum
      if (bcnums.find(bcnum) != bcnums.end()) {
        LOGF(info, "candCase %i bcnum %i allready found! ", candCase, bcnum);
        registry.get<TH1>(HIST("candCase"))->Fill(4, 1.);
        return;
      } else {
        bcnums.insert(bcnum);
      }

      // update histogram nDGperRun
      registry.get<TH1>(HIST("nDGperRun"))->Fill(Form("%d", run), 1.);
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

        registry.get<TH1>(HIST("TPCChi2NCl1"))->Fill(trk1.tpcChi2NCl(), 1.);
        registry.get<TH1>(HIST("TPCChi2NCl2"))->Fill(trk2.tpcChi2NCl(), 1.);
        registry.get<TH2>(HIST("2Trackpt1eta1"))->Fill(trk1.pt(), v1.Eta(), 1.);
        registry.get<TH2>(HIST("2Trackpt2eta2"))->Fill(trk2.pt(), v2.Eta(), 1.);
        registry.get<TH2>(HIST("2Trackpt1pt2"))->Fill(trk1.pt(), trk2.pt(), 1.);
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
        /*
        auto mom = sqrt(pow(track.px(), 2) + pow(track.py(), 2) + pow(track.pz(), 2));
        registry.get<TH2>(HIST("nSigmaTPCPEl"))->Fill(track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.get<TH2>(HIST("nSigmaTPCPPi"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.get<TH2>(HIST("nSigmaTPCPMu"))->Fill(track.tpcInnerParam(), track.tpcNSigmaMu());
        registry.get<TH2>(HIST("nSigmaTPCPKa"))->Fill(track.tpcInnerParam(), track.tpcNSigmaKa());
        registry.get<TH2>(HIST("nSigmaTPCPPr"))->Fill(track.tpcInnerParam(), track.tpcNSigmaPr());
        */
        registry.get<TH2>(HIST("nSigmaTPCPEl"))->Fill(track.pt(), track.tpcNSigmaEl());
        registry.get<TH2>(HIST("nSigmaTPCPPi"))->Fill(track.pt(), track.tpcNSigmaPi());
        registry.get<TH2>(HIST("nSigmaTPCPMu"))->Fill(track.pt(), track.tpcNSigmaMu());
        registry.get<TH2>(HIST("nSigmaTPCPKa"))->Fill(track.pt(), track.tpcNSigmaKa());
        registry.get<TH2>(HIST("nSigmaTPCPPr"))->Fill(track.pt(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          LOGF(debug, "tofNSigmaPi %f", track.tofNSigmaPi());
          /*
          registry.get<TH2>(HIST("nSigmaTOFPEl"))->Fill(mom, track.tofNSigmaEl());
          registry.get<TH2>(HIST("nSigmaTOFPPi"))->Fill(mom, track.tofNSigmaPi());
          registry.get<TH2>(HIST("nSigmaTOFPMu"))->Fill(mom, track.tofNSigmaMu());
          registry.get<TH2>(HIST("nSigmaTOFPKa"))->Fill(mom, track.tofNSigmaKa());
          registry.get<TH2>(HIST("nSigmaTOFPPr"))->Fill(mom, track.tofNSigmaPr());
          */
          registry.get<TH2>(HIST("nSigmaTOFPEl"))->Fill(track.pt(), track.tofNSigmaEl());
          registry.get<TH2>(HIST("nSigmaTOFPPi"))->Fill(track.pt(), track.tofNSigmaPi());
          registry.get<TH2>(HIST("nSigmaTOFPMu"))->Fill(track.pt(), track.tofNSigmaMu());
          registry.get<TH2>(HIST("nSigmaTOFPKa"))->Fill(track.pt(), track.tofNSigmaKa());
          registry.get<TH2>(HIST("nSigmaTOFPPr"))->Fill(track.pt(), track.tofNSigmaPr());
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
