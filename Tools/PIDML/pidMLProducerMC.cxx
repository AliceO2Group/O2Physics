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

/// \file pidMLProducerMc.cxx
/// \brief Produce PID ML skimmed data from MC files.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/PIDML/pidML.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct PidMlProducerMc {
  Produces<aod::PidTracksMcMl> pidTracksTableML;
  Produces<aod::PidTracksMc> pidTracksTable;

  Filter trackFilter = requireGlobalTrackInFilter();

  using BigTracksML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;
  using MyCollisionML = aod::Collisions::iterator;
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;
  using MyCollision = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>::iterator;

  HistogramRegistry registry{
    "registry",
    {{"hTPCSigvsPt", "TPC signal vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC signal", {HistType::kTH2F, {{500, 0., 10.}, {1000, 0., 600.}}}},
     {"hTOFBetavsPt", "TOF beta vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF beta", {HistType::kTH2F, {{500, 0., 10.}, {500, 0., 2.}}}},
     {"hTRDSigvsPt", "TRD signal vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TRD signal", {HistType::kTH2F, {{500, 0., 10.}, {2500, 0., 100.}}}}}};

  void processML(MyCollisionML const& collision, BigTracksML const& tracks, aod::McParticles const& mctracks)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = (uint8_t)mcParticle.isPhysicalPrimary();
      pidTracksTableML(track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                       track.tofSignal(), track.beta(),
                       track.p(), track.pt(), track.px(), track.py(), track.pz(),
                       track.sign(),
                       track.x(), track.y(), track.z(),
                       track.alpha(),
                       track.trackType(),
                       track.tpcNClsShared(),
                       track.dcaXY(), track.dcaZ(),
                       mcParticle.pdgCode(),
                       isPrimary);

      registry.fill(HIST("hTPCSigvsPt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("hTOFBetavsPt"), track.pt(), track.beta());
      registry.fill(HIST("hTRDSigvsPt"), track.pt(), track.trdSignal());
    }
  }
  PROCESS_SWITCH(PidMlProducerMc, processML, "Produce only ML MC essential data", true);

  void processAll(MyCollision const& collision, BigTracks const& tracks, aod::McParticles const& mctracks)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = (uint8_t)mcParticle.isPhysicalPrimary();
      pidTracksTable(collision.centRun2V0M(),
                     collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                     collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                     collision.multZNA(), collision.multZNC(),
                     collision.multTracklets(), collision.multTPC(),
                     track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                     track.trackEtaEmcal(), track.trackPhiEmcal(),
                     track.tofSignal(), track.beta(),
                     track.p(), track.pt(), track.px(), track.py(), track.pz(),
                     track.sign(),
                     track.x(), track.y(), track.z(),
                     track.alpha(),
                     track.trackType(),
                     track.tpcNClsShared(),
                     track.dcaXY(), track.dcaZ(),
                     track.tpcNSigmaEl(), track.tpcExpSigmaEl(), track.tpcExpSignalDiffEl(),
                     track.tofNSigmaEl(), track.tofExpSigmaEl(), track.tofExpSignalDiffEl(),
                     track.tpcNSigmaMu(), track.tpcExpSigmaMu(), track.tpcExpSignalDiffMu(),
                     track.tofNSigmaMu(), track.tofExpSigmaMu(), track.tofExpSignalDiffMu(),
                     track.tpcNSigmaPi(), track.tpcExpSigmaPi(), track.tpcExpSignalDiffPi(),
                     track.tofNSigmaPi(), track.tofExpSigmaPi(), track.tofExpSignalDiffPi(),
                     track.tpcNSigmaKa(), track.tpcExpSigmaKa(), track.tpcExpSignalDiffKa(),
                     track.tofNSigmaKa(), track.tofExpSigmaKa(), track.tofExpSignalDiffKa(),
                     track.tpcNSigmaPr(), track.tpcExpSigmaPr(), track.tpcExpSignalDiffPr(),
                     track.tofNSigmaPr(), track.tofExpSigmaPr(), track.tofExpSignalDiffPr(),
                     mcParticle.pdgCode(),
                     isPrimary);

      registry.fill(HIST("hTPCSigvsPt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("hTOFBetavsPt"), track.pt(), track.beta());
      registry.fill(HIST("hTRDSigvsPt"), track.pt(), track.trdSignal());
    }
  }
  PROCESS_SWITCH(PidMlProducerMc, processAll, "Produce all MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMlProducerMc>(cfgc)};
}
