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
// \brief Analyses UD tables (DGCandidates, DGTracks) of DG candidates produced with DGCandProducer
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  06.06.2022

#include "PWGUD/Core/DGPIDSelector.h"
#include "PWGUD/Core/UDGoodRunSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <set>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional printouts"};
  Configurable<int> candCaseSel{"CandCase", {}, "0: all Cands, 1: only ColCands,2: only BCCands"};
  Configurable<std::string> goodRunsFile{"goodRunsFile", {}, "json with list of good runs"};

  // ccdb
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // a pdg object
  TDatabasePDG* pdg = nullptr;

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
    {}};

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;
  // using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  // a function to fill 2Prong histograms
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
    auto signalTPC1 = tr1.tpcSignal();
    auto signalTOF1 = tr1.tofSignal() / 1.E4;

    auto tr2 = dgtracks.begin() + ivm.trkinds()[1];
    auto m2 = particleMass(pdg, pidsel.getAnaPars().PIDs()[1]);
    auto ene2 = sqrt(pow(tr2.px(), 2.) + pow(tr2.py(), 2.) + pow(tr2.pz(), 2.) + m2);
    auto lv2 = TLorentzVector(tr2.px(), tr2.py(), tr2.pz(), ene2);
    auto signalTPC2 = tr2.tpcSignal();
    auto signalTOF2 = tr2.tofSignal() / 1.E4;

    LOGF(debug, "TOF signals %f %f", signalTOF1, signalTOF2);

    registry.fill(HIST("2Prong/TPCsignal1"), tr1.tpcInnerParam(), signalTPC1);
    registry.fill(HIST("2Prong/TPCsignal2"), tr2.tpcInnerParam(), signalTPC2);
    registry.fill(HIST("2Prong/sig1VsSig2TPC"), signalTPC1, signalTPC2);
    registry.fill(HIST("2Prong/eta1Vseta2"), lv1.Eta(), lv2.Eta());

    if (tr1.hasTOF()) {
      registry.fill(HIST("2Prong/TOFsignal1"), lv1.P(), signalTOF1);
    }
    if (tr2.hasTOF()) {
      registry.fill(HIST("2Prong/TOFsignal2"), lv2.P(), signalTOF2);
    }
    if (tr1.hasTOF() && tr2.hasTOF()) {
      registry.fill(HIST("2Prong/sig1VsSig2TOF"), signalTOF1, signalTOF2);
    }
  }

  void init(InitContext& context)
  {
    // initalise ccdb
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // PDG
    pdg = TDatabasePDG::Instance();

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

    if (context.mOptions.get<bool>("processReco")) {
      registry.add("stat/candCaseAll", "Types of all DG candidates", {HistType::kTH1F, {{3, -0.5, 2.5}}});
      registry.add("stat/candCaseSel", "Types of all selectedDG candidates", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("stat/nDGperRun", "Number of DG collisions per run", {HistType::kTH1D, {{1, 0, 1}}});
      registry.add("stat/nPVtracks", "Number of PV tracks of analyzed collisions", {HistType::kTH1D, {{51, -0.5, 50.5}}});

      registry.add("tracks/nSigmaTPCPEl", "nSigma TPC for electrons", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
      registry.add("tracks/nSigmaTPCPPi", "nSigma TPC for pions", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
      registry.add("tracks/nSigmaTPCPMu", "nSigma TPC for muons", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
      registry.add("tracks/nSigmaTPCPKa", "nSigma TPC for kaons", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});
      registry.add("tracks/nSigmaTPCPPr", "nSigma TPC for protons", {HistType::kTH2F, {axispt, {100, -20.0, 20.0}}});

      const AxisSpec axisnsTOF{nsTOFAxis, "nSigma TOF axis for histograms"};
      registry.add("tracks/nSigmaTOFPEl", "nSigma TOF for electrons versus pT", {HistType::kTH2F, {axispt, axisnsTOF}});
      registry.add("tracks/nSigmaTOFPPi", "nSigma TOF for pions versus pT", {HistType::kTH2F, {axispt, axisnsTOF}});
      registry.add("tracks/nSigmaTOFPMu", "nSigma TOF for muons versus pT", {HistType::kTH2F, {axispt, axisnsTOF}});
      registry.add("tracks/nSigmaTOFPKa", "nSigma TOF for kaons versus pT", {HistType::kTH2F, {axispt, axisnsTOF}});
      registry.add("tracks/nSigmaTOFPPr", "nSigma TOF for protons versus pT", {HistType::kTH2F, {axispt, axisnsTOF}});

      registry.add("tracks/trackHits", "Track hits in various detectors", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("tracks/dcaXYDG", "dcaXY in DG candidates", {HistType::kTH1F, {{100, -0.2, 0.2}}});
      registry.add("tracks/dcaZDG", "dcaZ in DG candidates", {HistType::kTH1F, {{100, -0.5, 0.5}}});
      registry.add("tracks/TPCNCl", "Number of found TPC clusters", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("tracks/TPCChi2NCl", "TPC chi2 per cluster of tracks", {HistType::kTH1F, {{200, 0., 50.}}});
      registry.add("tracks/ptTrkdcaXYDG", "dcaXY versus track pT in DG candidates", {HistType::kTH2F, {axispt, {100, -0.2, 0.2}}});
      registry.add("tracks/ptTrkdcaZDG", "dcaZ versus track pT in DG candidates", {HistType::kTH2F, {axispt, {100, -0.5, 0.5}}});

      registry.add("system/nUnlikeIVMs", "Number of IVMs per DG collision", {HistType::kTH1F, {{36, -0.5, 35.5}}});
      registry.add("system/unlikeIVMptSysDG", "Invariant mass versus system pT in DG candidates", {HistType::kTH2F, {axisIVM, axispt}});
      registry.add("system/unlikeIVMptTrkDG", "Invariant mass versus track pT in DG candidates", {HistType::kTH2F, {axisIVM, axispt}});
      registry.add("system/nLikeIVMs", "Number of IVMs per DG collision", {HistType::kTH1F, {{36, -0.5, 35.5}}});
      registry.add("system/likeIVMptSysDG", "Invariant mass versus system pT in DG candidates", {HistType::kTH2F, {axisIVM, axispt}});
      registry.add("system/likeIVMptTrkDG", "Invariant mass versus track pT in DG candidates", {HistType::kTH2F, {axisIVM, axispt}});

      // FIT signals
      registry.add("FIT/FT0AAmplitude", "Total amplitude in FV0A", {HistType::kTH1F, {{5000, 0., 5000.}}});
      registry.add("FIT/FT0CAmplitude", "Total amplitude in FT0A", {HistType::kTH1F, {{5000, 0., 5000.}}});
      registry.add("FIT/FV0AAmplitude", "Total amplitude in FT0C", {HistType::kTH1F, {{5000, 0., 5000.}}});
      registry.add("FIT/FDDAAmplitude", "Total amplitude in FDDA", {HistType::kTH1F, {{5000, 0., 5000.}}});
      registry.add("FIT/FDDCAmplitude", "Total amplitude in FDDC", {HistType::kTH1F, {{5000, 0., 5000.}}});

      registry.add("FIT/BBFV0A", "FV0A signal in neighbouring BCs", {HistType::kTH1F, {{32, -16.5, 15.5}}});
      registry.add("FIT/BBFT0A", "FT0A signal in neighbouring BCs", {HistType::kTH1F, {{32, -16.5, 15.5}}});
      registry.add("FIT/BBFT0C", "FT0C signal in neighbouring BCs", {HistType::kTH1F, {{32, -16.5, 15.5}}});
      registry.add("FIT/BBFDDA", "FDDA signal in neighbouring BCs", {HistType::kTH1F, {{32, -16.5, 15.5}}});
      registry.add("FIT/BBFDDC", "FDDC signal in neighbouring BCs", {HistType::kTH1F, {{32, -16.5, 15.5}}});

      // 2 track events
      registry.add("2Prong/TPCNCl1", "Number of found TPC clusters of track 1", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("2Prong/TPCNCl2", "Number of found TPC clusters of track 2", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("2Prong/TPCChi2NCl1", "TPC chi2 of track 1", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("2Prong/TPCChi2NCl2", "TPC chi2 of track 2", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("2Prong/TPCsignal1", "TPC signal of track 1", {HistType::kTH2F, {{1000, 0., 10.}, {5000, 0., 500.}}});
      registry.add("2Prong/TPCsignal2", "TPC signal of track 2", {HistType::kTH2F, {{1000, 0., 10.}, {5000, 0., 500.}}});
      registry.add("2Prong/sig1VsSig2TPC", "TPC signals of track 1 versus track 2", {HistType::kTH2F, {{300, 0., 150.}, {300, 0., 150.}}});
      registry.add("2Prong/TOFsignal1", "TOF signal of track 1", {HistType::kTH2F, {{1000, 0., 10.}, {2000, -5., 5.}}});
      registry.add("2Prong/TOFsignal2", "TOF signal of track 2", {HistType::kTH2F, {{1000, 0., 10.}, {2000, -5., 5.}}});
      registry.add("2Prong/sig1VsSig2TOF", "TOF signals of track 1 versus track 2", {HistType::kTH2F, {{1000, -5., 5.}, {1000, -5., 5.}}});
      registry.add("2Prong/eta1Vseta2", "etas of track 1 versus track 2", {HistType::kTH2F, {{200, -2.0, 2.0}, {200, -2.0, 2.0}}});
      registry.add("2Prong/pt1pt2", "pTs of track 1 versus track 2", {HistType::kTH2F, {axispt, axispt}});
      registry.add("2Prong/pt1eta1", "pT versus eta of track 1", {HistType::kTH2F, {axispt, {200, -2.0, 2.0}}});
      registry.add("2Prong/pt2eta2", "pT versus eta of track 2", {HistType::kTH2F, {axispt, {200, -2.0, 2.0}}});
      registry.add("2Prong/Angle", "Angle between both tracks", {HistType::kTH1F, {{175, -0.2, 3.3}}});
      registry.add("2Prong/AngleIVM", "Angle versis invariant mass", {HistType::kTH2F, {axisIVM, {175, -0.2, 3.3}}});
      registry.add("2Prong/pt1IVM", "pT of track 1 versus invariant mass", {HistType::kTH2F, {axisIVM, axispt}});
      registry.add("2Prong/pt2IVM", "pT of track 2 versus invariant mass", {HistType::kTH2F, {axisIVM, axispt}});
      registry.add("2Prong/eta1IVM", "eta of track 1 versus invariant mass", {HistType::kTH2F, {axisIVM, {200, -2.0, 2.0}}});
      registry.add("2Prong/eta2IVM", "eta of track 2 versus invariant mass", {HistType::kTH2F, {axisIVM, {200, -2.0, 2.0}}});
      registry.add("2Prong/chi2NCl1IVM", "TPC chi2 of track 1 versus invariant mass", {HistType::kTH2F, {axisIVM, {100, 0, 5.0}}});
      registry.add("2Prong/chi2NCl2IVM", "TPC chi2 of track 2 versus invariant mass", {HistType::kTH2F, {axisIVM, {100, 0, 5.0}}});
      registry.add("2Prong/NCl1IVM", "Number of found TPC clusters of track 1 versus invariant mass", {HistType::kTH2F, {axisIVM, {200, 0, 200.}}});
      registry.add("2Prong/NCl2IVM", "Number of found TPC clusters of track 2 versus invariant mass", {HistType::kTH2F, {axisIVM, {200, 0, 200.}}});
    }

    if (context.mOptions.get<bool>("processMcTruth")) {
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {axisIVM, axispt}});
    }
  }

  // PV contributors
  Filter PVContributorFilter = aod::udtrack::isPVContributor == true;
  using PVTracks = soa::Filtered<UDTracksFull>;

  void processReco(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks, PVTracks const& PVContributors)
  {
    // count collisions
    registry.fill(HIST("stat/candCaseAll"), 0., 1.);
    registry.fill(HIST("stat/nPVtracks"), dgcand.numContrib(), 1.);

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
    registry.fill(HIST("stat/candCaseAll"), 1, 1.);

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
    } else if (dgcand.posX() == -3. && dgcand.posY() == 3. && dgcand.posZ() == -3.) {
      candCase = 3;
    }
    if (candCaseSel > 0 && candCase != candCaseSel) {
      return;
    }
    registry.fill(HIST("stat/candCaseAll"), 2, 1.);

    // fill FIT amplitude histograms
    registry.fill(HIST("FIT/FT0AAmplitude"), dgcand.totalFT0AmplitudeA(), 1.);
    registry.fill(HIST("FIT/FT0CAmplitude"), dgcand.totalFT0AmplitudeC(), 1.);
    registry.fill(HIST("FIT/FV0AAmplitude"), dgcand.totalFV0AmplitudeA(), 1.);
    registry.fill(HIST("FIT/FDDAAmplitude"), dgcand.totalFDDAmplitudeA(), 1.);
    registry.fill(HIST("FIT/FDDCAmplitude"), dgcand.totalFDDAmplitudeC(), 1.);

    // skip events with too few/many tracks
    if (dgcand.numContrib() != PVContributors.size()) {
      LOGF(info, "Missmatch of PVContributors %d != %d", dgcand.numContrib(), PVContributors.size());
    }
    if (dgcand.numContrib() < anaPars.minNTracks() || dgcand.numContrib() > anaPars.maxNTracks()) {
      LOGF(debug, "Rejected 1: %d not in range [%d, %d].", dgcand.numContrib(), anaPars.minNTracks(), anaPars.maxNTracks());
      return;
    }

    // skip events with out-of-range net charge
    auto netChargeValues = anaPars.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), dgcand.netCharge()) == netChargeValues.end()) {
      LOGF(debug, "Rejected 2: %d not in set.", dgcand.netCharge());
      return;
    }

    // skip events with out-of-range rgtrwTOF
    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    auto minRgtrwTOF = candCase != 1 ? 1.0 : anaPars.minRgtrwTOF();
    if (rtrwTOF < minRgtrwTOF) {
      LOGF(debug, "Rejected 3: %f below threshold of %f.", rtrwTOF, minRgtrwTOF);
      return;
    }

    // check FIT information
    auto bitMin = anaPars.dBCMin() + 16;
    auto bitMax = anaPars.dBCMax() + 16;
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (anaPars.FITvetoes()[0] && TESTBIT(dgcand.bbFV0Apf(), bit))
        return;
      if (anaPars.FITvetoes()[1] && TESTBIT(dgcand.bbFT0Apf(), bit))
        return;
      if (anaPars.FITvetoes()[2] && TESTBIT(dgcand.bbFT0Cpf(), bit))
        return;
      if (anaPars.FITvetoes()[3] && TESTBIT(dgcand.bbFDDApf(), bit))
        return;
      if (anaPars.FITvetoes()[4] && TESTBIT(dgcand.bbFDDCpf(), bit))
        return;
    }

    // fill BBFlag histograms
    for (auto bit = 0; bit < 33; bit++) {
      registry.fill(HIST("FIT/BBFV0A"), bit - 16, TESTBIT(dgcand.bbFV0Apf(), bit));
      registry.fill(HIST("FIT/BBFT0A"), bit - 16, TESTBIT(dgcand.bbFT0Apf(), bit));
      registry.fill(HIST("FIT/BBFT0C"), bit - 16, TESTBIT(dgcand.bbFT0Cpf(), bit));
      registry.fill(HIST("FIT/BBFDDA"), bit - 16, TESTBIT(dgcand.bbFDDApf(), bit));
      registry.fill(HIST("FIT/BBFDDC"), bit - 16, TESTBIT(dgcand.bbFDDCpf(), bit));
    }

    // find track combinations which are compatible with PID cuts
    auto nIVMs = pidsel.computeIVMs(PVContributors);

    // process the unlike sign combinations
    if (nIVMs[0] == 0 && nIVMs[1] == 0) {
      LOGF(debug, "Rejected 4: no IVMs.");
      return;
    }
    LOGF(debug, "nIVMs %d / %d", nIVMs[0], nIVMs[1]);

    // update histogram stat/candCase and stat/nDGperRun
    registry.fill(HIST("stat/candCaseSel"), 0, 1.);
    registry.fill(HIST("stat/candCaseSel"), candCase, 1.);

    /*
    // check bcnum
    if (bcnums.find(bcnum) != bcnums.end()) {
      LOGF(debug, "candCase %d bcnum %d allready found! ", candCase, bcnum);
      registry.fill(HIST("stat/candCaseSel"), 4, 1.);
    } else {
      bcnums.insert(bcnum);
    }
    */
    registry.get<TH1>(HIST("stat/nDGperRun"))->Fill(Form("%d", run), 1);

    // update histograms
    int goodIVMs = 0;
    for (auto ivm : pidsel.unlikeIVMs()) {
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
        }

        // update 2Prong histograms
        registry.fill(HIST("2Prong/Angle"), angle, 1.);
        registry.fill(HIST("2Prong/AngleIVM"), ivm.M(), angle, 1.);

        registry.fill(HIST("2Prong/TPCChi2NCl1"), trk1.tpcChi2NCl(), 1.);
        registry.fill(HIST("2Prong/TPCChi2NCl2"), trk2.tpcChi2NCl(), 1.);
        registry.fill(HIST("2Prong/pt1eta1"), trk1.pt(), v1.Eta(), 1.);
        registry.fill(HIST("2Prong/pt2eta2"), trk2.pt(), v2.Eta(), 1.);
        registry.fill(HIST("2Prong/pt1pt2"), trk1.pt(), trk2.pt(), 1.);
        registry.fill(HIST("2Prong/pt1IVM"), ivm.M(), trk1.pt(), 1.);
        registry.fill(HIST("2Prong/pt2IVM"), ivm.M(), trk2.pt(), 1.);
        registry.fill(HIST("2Prong/eta1IVM"), ivm.M(), v1.Eta(), 1.);
        registry.fill(HIST("2Prong/eta2IVM"), ivm.M(), v2.Eta(), 1.);
        registry.fill(HIST("2Prong/chi2NCl1IVM"), ivm.M(), trk1.tpcChi2NCl(), 1.);
        registry.fill(HIST("2Prong/chi2NCl2IVM"), ivm.M(), trk2.tpcChi2NCl(), 1.);

        auto nTPCCL = trk1.tpcNClsFindable() - trk1.tpcNClsFindableMinusFound();
        registry.fill(HIST("2Prong/TPCNCl1"), nTPCCL, 1.);
        registry.fill(HIST("2Prong/NCl1IVM"), ivm.M(), nTPCCL, 1.);
        nTPCCL = trk2.tpcNClsFindable() - trk2.tpcNClsFindableMinusFound();
        registry.fill(HIST("2Prong/TPCNCl2"), nTPCCL, 1.);
        registry.fill(HIST("2Prong/NCl2IVM"), ivm.M(), nTPCCL, 1.);

        fillSignalHists(ivm, PVContributors, pidsel);
      }
      goodIVMs++;

      // update system/IVMptSysDG
      registry.fill(HIST("system/unlikeIVMptSysDG"), ivm.M(), ivm.Perp());

      // loop over tracks of IVM and fill related histograms
      for (auto ind : ivm.trkinds()) {
        auto track = PVContributors.begin() + ind;
        registry.fill(HIST("system/unlikeIVMptTrkDG"), ivm.M(), track.pt());

        registry.fill(HIST("tracks/trackHits"), 0., 1.);
        registry.fill(HIST("tracks/trackHits"), 1., track.hasITS() * 1.);
        registry.fill(HIST("tracks/trackHits"), 2., track.hasTPC() * 1.);
        registry.fill(HIST("tracks/trackHits"), 3., track.hasTRD() * 1.);
        registry.fill(HIST("tracks/trackHits"), 4., track.hasTOF() * 1.);

        registry.fill(HIST("tracks/dcaXYDG"), track.dcaXY());
        registry.fill(HIST("tracks/ptTrkdcaXYDG"), track.pt(), track.dcaXY());
        registry.fill(HIST("tracks/dcaZDG"), track.dcaZ());
        registry.fill(HIST("tracks/ptTrkdcaZDG"), track.pt(), track.dcaZ());
        registry.fill(HIST("tracks/TPCNCl"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound(), 1.);
        registry.fill(HIST("tracks/TPCChi2NCl"), track.tpcChi2NCl(), 1.);

        // fill nSigma histograms
        /*
        auto mom = sqrt(pow(track.px(), 2) + pow(track.py(), 2) + pow(track.pz(), 2));
        registry.fill(HIST("tracks/nSigmaTPCPEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.fill(HIST("tracks/nSigmaTPCPPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.fill(HIST("tracks/nSigmaTPCPMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
        registry.fill(HIST("tracks/nSigmaTPCPKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        registry.fill(HIST("tracks/nSigmaTPCPPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        */
        registry.fill(HIST("tracks/nSigmaTPCPEl"), track.pt(), track.tpcNSigmaEl());
        registry.fill(HIST("tracks/nSigmaTPCPPi"), track.pt(), track.tpcNSigmaPi());
        registry.fill(HIST("tracks/nSigmaTPCPMu"), track.pt(), track.tpcNSigmaMu());
        registry.fill(HIST("tracks/nSigmaTPCPKa"), track.pt(), track.tpcNSigmaKa());
        registry.fill(HIST("tracks/nSigmaTPCPPr"), track.pt(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          LOGF(debug, "tofNSigmaPi %f", track.tofNSigmaPi());
          /*
          registry.fill(HIST("nSigmaTOFPEl"), mom, track.tofNSigmaEl());
          registry.fill(HIST("nSigmaTOFPPi"), mom, track.tofNSigmaPi());
          registry.fill(HIST("nSigmaTOFPMu"), mom, track.tofNSigmaMu());
          registry.fill(HIST("nSigmaTOFPKa"), mom, track.tofNSigmaKa());
          registry.fill(HIST("nSigmaTOFPPr"), mom, track.tofNSigmaPr());
          */
          registry.fill(HIST("tracks/nSigmaTOFPEl"), track.pt(), track.tofNSigmaEl());
          registry.fill(HIST("tracks/nSigmaTOFPPi"), track.pt(), track.tofNSigmaPi());
          registry.fill(HIST("tracks/nSigmaTOFPMu"), track.pt(), track.tofNSigmaMu());
          registry.fill(HIST("tracks/nSigmaTOFPKa"), track.pt(), track.tofNSigmaKa());
          registry.fill(HIST("tracks/nSigmaTOFPPr"), track.pt(), track.tofNSigmaPr());
        }
      }
    }
    LOGF(debug, "goodIVMs %d", goodIVMs);
    registry.fill(HIST("system/nUnlikeIVMs"), goodIVMs, 1.);

    // process the like sign combinations
    goodIVMs = 0;
    for (auto ivm : pidsel.likeIVMs()) {
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
        }
      }
      goodIVMs++;

      // update system/IVMptSysDG
      registry.fill(HIST("system/likeIVMptSysDG"), ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = PVContributors.begin() + ind;
        registry.fill(HIST("system/likeIVMptTrkDG"), ivm.M(), track.pt());
      }
    }
    registry.fill(HIST("system/nLikeIVMs"), goodIVMs, 1.);
  }

  PROCESS_SWITCH(DGCandAnalyzer, processReco, "Analyse reconstructed data", true);

  using UDMcCollisionsFull = aod::UDMcCollisions;
  using UDMcCollisionFull = UDMcCollisionsFull::iterator;
  using UDMcTracksFull = aod::UDMcParticles;

  void processMcTruth(UDMcCollisionFull const& /*mcCollison*/, UDMcTracksFull const& mcParts)
  {

    // which type of event is this
    // 0: MB
    // 1: Pythia diffractive
    // 2: GRANIITTI diffractive
    bool isPythiaDiff = udhelpers::isPythiaCDE(mcParts);
    bool isGraniittiDiff = udhelpers::isGraniittiCDE(mcParts);
    registry.get<TH1>(HIST("mcTruth/collType"))->Fill(0., 1.);
    registry.get<TH1>(HIST("mcTruth/collType"))->Fill(1., (!isPythiaDiff && !isGraniittiDiff) * 1.);
    registry.get<TH1>(HIST("mcTruth/collType"))->Fill(2., isPythiaDiff * 1.);
    registry.get<TH1>(HIST("mcTruth/collType"))->Fill(3., isGraniittiDiff * 1.);

    // compute GRANIITTI event invariant mass
    if (!isGraniittiDiff) {
      return;
    }
    auto ivm = udhelpers::ivmGraniittiCDE(mcParts);

    // update histograms
    registry.get<TH2>(HIST("mcTruth/IVMpt"))->Fill(ivm.M(), ivm.Perp(), 1.);
  }

  PROCESS_SWITCH(DGCandAnalyzer, processMcTruth, "Analyse MC truth", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandAnalyzer>(cfgc, TaskName{"dgcandanalyzer"}),
  };
}
