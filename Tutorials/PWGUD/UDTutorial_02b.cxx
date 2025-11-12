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

#include "PWGUD/Core/DGPIDSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDTutorial02b {

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  ConfigurableAxis IVMAxis{"IVMAxis", {350, 0.0, 3.5}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  ConfigurableAxis nsTPCAxis{"nsTPCAxis", {100, -20.0, 20.0}, "nSigma TPC axis"};
  ConfigurableAxis nsTOFAxis{"nsTOFAxis", {100, -100.0, 100.0}, "nSigma TOF axis"};

  // analysis cuts
  DGAnaparHolder anaPars = DGAnaparHolder();
  Configurable<DGAnaparHolder> DGPars{"anaPars", {}, "Analysis parameters"};

  // PID selector
  DGPIDSelector pidsel = DGPIDSelector();

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    anaPars = (DGAnaparHolder)DGPars;
    pidsel.init(anaPars);
    if (verbose) {
      LOGF(info, "<UDTutorial02>");
      pidsel.Print();
    }

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
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::TracksDCA, aod::UDTracksFlags>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {

    // skip events with too few/many tracks
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);
    if (dgcand.numContrib() < anaPars.minNTracks() || dgcand.numContrib() > anaPars.maxNTracks()) {
      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate rejected: Number of PV contributors %d not in range [%d, %d].", dgcand.numContrib(), anaPars.minNTracks(), anaPars.maxNTracks());
      }
      return;
    }

    // skip events with out-of-range net charge
    auto netChargeValues = anaPars.unlikeCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), dgcand.netCharge()) == netChargeValues.end()) {
      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate rejected: Net charge %d not in selected set.", dgcand.netCharge());
      }
      return;
    }

    // skip events with out-of-range rgtrwTOF (fraction-of-good-tracks-with-TOF-hit)
    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    if (rtrwTOF < anaPars.minRgtrwTOF()) {
      if (verbose) {
        LOGF(debug, "<UDTutorials02> Candidate rejected: rtrwTOF %f below threshold of %f.", rtrwTOF, anaPars.minRgtrwTOF());
      }
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

    // find track combinations which are compatible with PID cuts
    auto nIVMs = pidsel.computeIVMs(PVContributors);
    if (verbose) {
      LOGF(info, "<<UDTutorials02> Found %d/%d candidates", nIVMs[0], nIVMs[1]);
    }

    // update histograms
    for (auto ivm : pidsel.unlikeIVMs()) {
      // cut on pt-system
      if (ivm.Perp() < anaPars.minptsys() || ivm.Perp() > anaPars.maxptsys()) {
        continue;
      }

      if (verbose) {
        LOGF(info, "<UDTutorials02> Candidate accepted!");
      }
      registry.get<TH2>(HIST("dgcandidates/IVMptsys"))->Fill(ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = PVContributors.begin() + ind;
        registry.get<TH2>(HIST("dgcandidates/IVMpttrk"))->Fill(ivm.M(), track.pt());

        // fill nSigma histograms
        auto mom = sqrt(pow(track.px(), 2) + pow(track.py(), 2) + pow(track.pz(), 2));
        registry.get<TH2>(HIST("dgcandidates/nSigmaTPCEl"))->Fill(mom, track.tpcNSigmaEl());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTPCPi"))->Fill(mom, track.tpcNSigmaPi());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTPCMu"))->Fill(mom, track.tpcNSigmaMu());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTPCKa"))->Fill(mom, track.tpcNSigmaKa());
        registry.get<TH2>(HIST("dgcandidates/nSigmaTPCPr"))->Fill(mom, track.tpcNSigmaPr());
        if (track.hasTOF()) {
          LOGF(debug, "tofNSigmaPi %f", track.tofNSigmaPi());
          registry.get<TH2>(HIST("dgcandidates/nSigmaTOFEl"))->Fill(mom, track.tofNSigmaEl());
          registry.get<TH2>(HIST("dgcandidates/nSigmaTOFPi"))->Fill(mom, track.tofNSigmaPi());
          registry.get<TH2>(HIST("dgcandidates/nSigmaTOFMu"))->Fill(mom, track.tofNSigmaMu());
          registry.get<TH2>(HIST("dgcandidates/nSigmaTOFKa"))->Fill(mom, track.tofNSigmaKa());
          registry.get<TH2>(HIST("dgcandidates/nSigmaTOFPr"))->Fill(mom, track.tofNSigmaPr());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial02b>(cfgc, TaskName{"udtutorial02b"}),
  };
}
