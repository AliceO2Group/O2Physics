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
// \brief UD tutorial
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  April 2023

#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/PIDResponseTOF.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TVector3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDTutorial01 {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  ConfigurableAxis etaAxis{"etaAxis", {300, -1.5, 1.5}, ""};
  ConfigurableAxis sigTPCAxis{"sigTPCAxis", {100, -100.0, 100.0}, ""};
  ConfigurableAxis sigTOFAxis{"sigTOFAxis", {100, -100.0, 100.0}, ""};

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // Collision histograms
    registry.add("collisions/BC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    registry.add("collisions/multiplicityAll", "Multiplicity of all tracks; Tracks; Tracks", {HistType::kTH1F, {{201, -0.5, 200.5}}});
    registry.add("collisions/multiplicityPVC", "Multiplicity of PV contributors; PV contributors; Tracks", {HistType::kTH1F, {{201, -0.5, 200.5}}});

    // track histograms
    const AxisSpec axispt{ptAxis, "p_{T} axis"};
    const AxisSpec axiseta{etaAxis, "pseudo rapidity axis"};
    registry.add("tracks/QCAll", "Track QC of all tracks; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/QCPVC", "Track QC of PV contributors; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/ptAll", "track pt of all tracks; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axispt}});
    registry.add("tracks/ptPVC", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axispt}});
    registry.add("tracks/etavsptAll", "track eta versus pt of all tracks; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});
    registry.add("tracks/etavsptPVC", "track eta versus pt of PV contributors; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});

    const AxisSpec axisp{ptAxis, "momentum axis"};
    const AxisSpec axisTPCsig{sigTPCAxis, "TPC signal"};
    const AxisSpec axisTOFsig{sigTOFAxis, "TOF signal"};
    registry.add("tracks/TPCSignalvspAll", "TPC signal versus track momentum of all tracks; Track momentum [GeV/c]; TPC signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTPCsig}});
    registry.add("tracks/TPCSignalvspPVC", "TPC signal versus track momentum of PV contributors; Track momentum [GeV/c]; TPC signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTPCsig}});
    registry.add("tracks/TOFSignalvspAll", "TOF signal versus track momentum of all tracks; Track momentum [GeV/c]; TOF signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTOFsig}});
    registry.add("tracks/TOFSignalvspPVC", "TOF signal versus track momentum of PV contributors; Track momentum [GeV/c]; TOF signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTOFsig}});

    // FIT histograms
    registry.add("FIT/BBFV0A", "Beam-beam in V0A; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFT0A", "Beam-beam in T0A; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFT0C", "Beam-beam in T0C; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFDDA", "Beam-beam in FDDA; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFDDC", "Beam-beam in FDDA; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
  }

  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    if (verbose) {
      LOGF(info, "<UDTutorial01> DG candidate %d", dgcand.globalIndex());
    }

    // fill collision histograms
    registry.get<TH1>(HIST("collisions/multiplicityAll"))->Fill(dgtracks.size(), 1.);
    // select PV contributors
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);
    registry.get<TH1>(HIST("collisions/multiplicityPVC"))->Fill(PVContributors.size(), 1.);

    // relative BC number
    auto bcnum = dgcand.globalBC() % o2::constants::lhc::LHCMaxBunches;
    registry.get<TH1>(HIST("collisions/BC"))->Fill(bcnum, 1.);

    // fill track histograms
    if (verbose) {
      LOGF(info, "<UDTutorial01>   Number of tracks %d", dgtracks.size());
      LOGF(info, "<UDTutorial01>   Number of PV contributors %d", PVContributors.size());
    }
    for (auto track : dgtracks) {
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(0., 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(1., track.hasITS() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(2., track.hasTPC() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(3., track.hasTRD() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(4., track.hasTOF() * 1.);

      auto vtrk = TVector3(track.px(), track.py(), track.pz());
      registry.get<TH1>(HIST("tracks/ptAll"))->Fill(track.pt(), 1.);
      registry.get<TH2>(HIST("tracks/etavsptAll"))->Fill(vtrk.Eta(), track.pt(), 1.);

      auto signalTPC = track.tpcSignal() * track.sign();
      registry.get<TH2>(HIST("tracks/TPCSignalvspAll"))->Fill(vtrk.Mag(), signalTPC, 1.);
      auto signalTOF = track.tofSignal() * track.sign() / 1.E3;
      registry.get<TH2>(HIST("tracks/TOFSignalvspAll"))->Fill(vtrk.Mag(), signalTOF, 1.);

      if (track.isPVContributor()) {
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(0., 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(1., track.hasITS() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(2., track.hasTPC() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(3., track.hasTRD() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(4., track.hasTOF() * 1.);
        registry.get<TH1>(HIST("tracks/ptPVC"))->Fill(track.pt(), 1.);
        registry.get<TH2>(HIST("tracks/etavsptPVC"))->Fill(vtrk.Eta(), track.pt(), 1.);
        registry.get<TH2>(HIST("tracks/TPCSignalvspPVC"))->Fill(vtrk.Mag(), signalTPC, 1.);
        registry.get<TH2>(HIST("tracks/TOFSignalvspPVC"))->Fill(vtrk.Mag(), signalTOF, 1.);
      }
    }

    // fill FIT histograms
    for (auto bit = 0; bit < 33; bit++) {
      registry.get<TH1>(HIST("FIT/BBFV0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFV0Apf(), bit));
      registry.get<TH1>(HIST("FIT/BBFT0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Apf(), bit));
      registry.get<TH1>(HIST("FIT/BBFT0C"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Cpf(), bit));
      registry.get<TH1>(HIST("FIT/BBFDDA"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDApf(), bit));
      registry.get<TH1>(HIST("FIT/BBFDDC"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDCpf(), bit));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial01>(cfgc, TaskName{"udtutorial01"}),
  };
}
