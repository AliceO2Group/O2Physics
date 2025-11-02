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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDTutorial02a {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  ConfigurableAxis IVMAxis{"IVMAxis", {350, 0.0, 3.5}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  ConfigurableAxis nsTPCAxis{"nsTPCAxis", {100, -20.0, 20.0}, "nSigma TPC axis"};
  ConfigurableAxis nsTOFAxis{"nsTOFAxis", {100, -100.0, 100.0}, "nSigma TOF axis"};

  // a pdg object
  TDatabasePDG* pdg = nullptr;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    pdg = TDatabasePDG::Instance();

    // dgcandidates histograms
    const AxisSpec axisIVM{IVMAxis, "IVM axis"};
    const AxisSpec axispt{ptAxis, "pt axis"};
    registry.add("dgcandidates/IVMptsys", "IVM versus system pT; Invariant mass [GeV/c^{2}]; p_{T, system} {GeV/c]", {HistType::kTH2F, {axisIVM, axispt}});
    registry.add("dgcandidates/IVMpttrk", "IVM versus track pT; Invariant mass [GeV/c^{2}]; p_{T, track} {GeV/c]", {HistType::kTH2F, {axisIVM, axispt}});

    const AxisSpec axisnsTPC{nsTOFAxis, "nSigma TPC axis"};
    registry.add("dgcandidates/nSigmaTPCEl", "TPC nSigma electrons versus track pT; Track p_{T} [GeV/c]; TPC nSigma_{e}", {HistType::kTH2F, {axispt, axisnsTPC}});
    registry.add("dgcandidates/nSigmaTPCPi", "TPC nSigma pions versus track pT; Track p_{T} [GeV/c]; TPC nSigma_{pi}", {HistType::kTH2F, {axispt, axisnsTPC}});
    registry.add("dgcandidates/nSigmaTPCMu", "TPC nSigma muons versus track pT; Track p_{T} [GeV/c]; TPC nSigma_{mu}", {HistType::kTH2F, {axispt, axisnsTPC}});
    registry.add("dgcandidates/nSigmaTPCKa", "TPC nSigma kaon versus track pT; Track p_{T} [GeV/c]; TPC nSigma_{K}", {HistType::kTH2F, {axispt, axisnsTPC}});
    registry.add("dgcandidates/nSigmaTPCPr", "TPC nSigma protons versus track pT; Track p_{T} [GeV/c]; TPC nSigma_{p}", {HistType::kTH2F, {axispt, axisnsTPC}});

    const AxisSpec axisnsTOF{nsTOFAxis, "nSigma TOF axis"};
    registry.add("dgcandidates/nSigmaTOFEl", "TOF nSigma electrons versus track pT; Track p_{T} [GeV/c]; TOF nSigma_{e}", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("dgcandidates/nSigmaTOFPi", "TOF nSigma pions versus track pT; Track p_{T} [GeV/c]; TOF nSigma_{pi}", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("dgcandidates/nSigmaTOFMu", "TOF nSigma muons versus track pT; Track p_{T} [GeV/c]; TOF nSigma_{mu}", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("dgcandidates/nSigmaTOFKa", "TOF nSigma kaons versus track pT; Track p_{T} [GeV/c]; TOF nSigma_{K}", {HistType::kTH2F, {axispt, axisnsTOF}});
    registry.add("dgcandidates/nSigmaTOFPr", "TOF nSigma protons versus track pT; Track p_{T} [GeV/c]; TOF nSigma_{p}", {HistType::kTH2F, {axispt, axisnsTOF}});
  }

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {

    // skip events with too few/many tracks
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);
    if (dgcand.numContrib() != 2) {
      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate rejected: Number of PV contributors is %d", dgcand.numContrib());
      }
      return;
    }

    // skip events with net charge != 0
    if (dgcand.netCharge() != 0) {
      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate rejected: Net charge is %d", dgcand.netCharge());
      }
      return;
    }

    // skip events with out-of-range rgtrwTOF (fraction-of-good-tracks-with-TOF-hit)
    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    if (rtrwTOF < 0.5) {
      if (verbose) {
        LOGF(debug, "<UDTutorials02> Candidate rejected: rtrwTOF is %f", rtrwTOF);
      }
      return;
    }

    // check FIT information
    auto bitMin = -1 + 16;
    auto bitMax = 1 + 16;
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (TESTBIT(dgcand.bbFT0Apf(), bit) ||
          TESTBIT(dgcand.bbFT0Cpf(), bit) ||
          TESTBIT(dgcand.bbFV0Apf(), bit) ||
          TESTBIT(dgcand.bbFDDApf(), bit) ||
          TESTBIT(dgcand.bbFDDCpf(), bit)) {
        return;
      }
    }

    // check PID of tracks, use nSigmaTPC
    // cut on track pT
    for (auto trk : PVContributors) {
      if (trk.tpcNSigmaPi() < -3. || trk.tpcNSigmaPi() > 3.) {
        if (verbose) {
          LOGF(info, "<UDTutorials02> Candidate rejected: nSigmaTPC pion is %f", trk.tpcNSigmaPi());
        }
        return;
      }
      if (trk.pt() < 0.1) {
        if (verbose) {
          LOGF(info, "<UDTutorials02> Candidate rejected: Track pT is %f", trk.pt());
        }
        return;
      }
    }

    // compute invariant mass
    TParticlePDG* pion = pdg->GetParticle(211);
    TLorentzVector lvtmp;
    auto ivm = TLorentzVector(0., 0., 0., 0.);
    for (auto trk : PVContributors) {
      lvtmp.SetXYZM(trk.px(), trk.py(), trk.pz(), pion->Mass());
      ivm += lvtmp;
    }

    // cut on system pT
    if (ivm.Perp() < 0.1) {
      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate rejected: System pT is %f", ivm.Perp());
      }
      return;
    }

    // update histograms
    if (verbose) {
      LOGF(info, "<UDTutorials02> Candidate accepted!");
    }
    registry.get<TH2>(HIST("dgcandidates/IVMptsys"))->Fill(ivm.M(), ivm.Perp());
    for (auto trk : PVContributors) {
      registry.get<TH2>(HIST("dgcandidates/IVMpttrk"))->Fill(ivm.M(), trk.pt());

      // fill nSigma histograms
      auto mom = sqrt(pow(trk.px(), 2) + pow(trk.py(), 2) + pow(trk.pz(), 2));
      registry.get<TH2>(HIST("dgcandidates/nSigmaTPCEl"))->Fill(mom, trk.tpcNSigmaEl());
      registry.get<TH2>(HIST("dgcandidates/nSigmaTPCPi"))->Fill(mom, trk.tpcNSigmaPi());
      registry.get<TH2>(HIST("dgcandidates/nSigmaTPCMu"))->Fill(mom, trk.tpcNSigmaMu());
      registry.get<TH2>(HIST("dgcandidates/nSigmaTPCKa"))->Fill(mom, trk.tpcNSigmaKa());
      registry.get<TH2>(HIST("dgcandidates/nSigmaTPCPr"))->Fill(mom, trk.tpcNSigmaPr());
      if (trk.hasTOF()) {
        registry.get<TH2>(HIST("dgcandidates/nSigmaTOFEl"))->Fill(mom, trk.tofNSigmaEl());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTOFPi"))->Fill(mom, trk.tofNSigmaPi());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTOFMu"))->Fill(mom, trk.tofNSigmaMu());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTOFKa"))->Fill(mom, trk.tofNSigmaKa());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTOFPr"))->Fill(mom, trk.tofNSigmaPr());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial02a>(cfgc, TaskName{"udtutorial02a"}),
  };
}
