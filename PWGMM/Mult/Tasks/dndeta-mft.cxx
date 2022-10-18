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

#include <cmath>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"
// for CCDB access
#include <CCDB/BasicCCDBManager.h>
#include <chrono>

#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec PtAxis = {1001, -0.005, 10.005};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec MultAxis = {301, -0.5, 300.5};
AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec PhiAxis = {600, 0, 2 * M_PI};
AxisSpec EtaAxis = {18, -4.6, -1.};

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct PseudorapidityDensityMFT {
  Service<TDatabasePDG> pdg;

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};

  Configurable<bool> useZDiffCut{"useZDiffCut", true, "use Z difference cut"};
  Configurable<float> maxZDiff{"maxZDiff", 1.0f, "max allowed Z difference for reconstruced collisions (cm)"};
  //-----need access to CCDB to get the reweighting histogram
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> path{"ccdb-path", "Users/s/sherrman/My/Object", "base path to the ccdb object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // the histogram has been previously stored in the CCDB
  TH1D* histoReweight = nullptr;
  std::vector<long long unsigned> ambTrackIds;
  int counter = 0;
  //------

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {{301, -0.5, 300.5}, ZAxis}}}, //
      {"TracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},                //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}},                         //
      {"TracksPhiZvtx", "; #varphi; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {PhiAxis, ZAxis}}},             //
      {"TracksPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}},                            //
      {"EventSelection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}}                                    //
    }                                                                                                             //
  };

  void init(InitContext&)
  {

    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Selected INEL>0");
    x->SetBinLabel(4, "Rejected");
    x->SetBinLabel(5, "Good BCs");
    x->SetBinLabel(6, "BCs with collisions");
    x->SetBinLabel(7, "BCs with pile-up/splitting");

    if (doprocessGen) {
      registry.add({"EventsNtrkZvtxGen", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {{301, -0.5, 300.5}, ZAxis}}});
      registry.add({"EventsNtrkZvtxGen_t", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {{301, -0.5, 300.5}, ZAxis}}});
      registry.add({"TracksEtaZvtxGen", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksEtaZvtxGen_t", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"TracksPhiZvtxGen", "; #varphi; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {PhiAxis, ZAxis}}}); //
      registry.add({"TracksToPartPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});             //
      registry.add({"TracksPtEtaGen", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"TracksPtEtaGen_t", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"EventEfficiency", "; status; events", {HistType::kTH1F, {{5, 0.5, 5.5}}}});
      registry.add({"NotFoundEventZvtx", " ; #it{z}_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}});
      registry.add({"EventsZposDiff", " ; Z_{rec} - Z_{gen} (cm)", {HistType::kTH1F, {DeltaZAxis}}});
      registry.add({"EventsSplitMult", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});

      auto heff = registry.get<TH1>(HIST("EventEfficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }
    if (doprocessMultReweight) {
      ccdb->setURL(url.value);
      // Enabling object caching, otherwise each call goes to the CCDB server
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      // Not later than now, will be replaced by the value of the train creation
      // This avoids that users can replace objects **while** a train is running
      ccdb->setCreatedNotAfter(nolaterthan.value);
      LOGF(info, "Getting object %s", path.value.data());
      histoReweight = ccdb->getForTimeStamp<TH1D>(path.value, nolaterthan.value);
      if (!histoReweight) {
        LOGF(fatal, "object not found!");
      }
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {

    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (useEvSel && ((bc.selection()[evsel::kIsBBT0A] & bc.selection()[evsel::kIsBBT0C]) != 0))) {
        registry.fill(HIST("EventSelection"), 5.);
        cols.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("EventSelection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("EventSelection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging, "Collect event sample stats", true);

  Partition<aod::MFTTracks> sample = (aod::fwdtrack::eta < -2.8f) && (aod::fwdtrack::eta > -3.2f);
  expressions::Filter atrackFilter = (aod::fwdtrack::bestCollisionId >= 0) &&
                                     (aod::fwdtrack::etas < -2.0f) &&
                                     (aod::fwdtrack::etas > -3.9f) &&
                                     (nabs(aod::fwdtrack::bestDCAXY) <= 2.f);

  void processAmbi(aod::AmbiguousMFTTracks const& atracks, soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::MFTTracks const& tracks)
  {
    counter++;
    ambTrackIds.clear();
    for (auto& track : atracks) {
      auto mfttrack = track.mfttrack();
      if (mfttrack.has_collision()) {
        ambTrackIds.push_back(track.mfttrackId());
      }
    }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processAmbi, "Fills vector of ambiTracks", true);

  void processMult(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks, soa::SmallGroups<soa::Join<aod::AmbiguousMFTTracks, aod::BestCollisionsFwd>> const& atracks)
  {

    if (tracks.size() == 0) {
      return;
    }
    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sample->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex());
      auto Ntrk = perCollisionSample.size() + atracks.size();

      registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);
      for (auto& track : atracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.etas(), z);
        float phi = track.phis();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksPhiEta"), phi, track.etas());
        registry.fill(HIST("TracksPhiZvtx"), phi, z);
        registry.fill(HIST("TracksPtEta"), track.pts(), track.etas());
      }

      for (auto& track : tracks) {
        if (find(ambTrackIds.begin(), ambTrackIds.end(), track.globalIndex()) != ambTrackIds.end()) {
          continue; // this track has already been reassigned to bestcollision, don't double count
        }
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksPhiEta"), phi, track.eta());
        registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
        if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
          registry.fill(HIST("TracksPhiZvtx"), phi, z);
        }
      }

    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMult, "Process reco or data info", true);

  void processMultReweight(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks, soa::SmallGroups<soa::Join<aod::AmbiguousMFTTracks, aod::BestCollisionsFwd>> const& atracks)
  {

    if (!doprocessGen) {
      LOGP(debug, "You can't enable processMultReweight if not analysing MC");
      return;
    }
    if (tracks.size() == 0) {
      return;
    }
    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sample->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex());
      auto Ntrk = perCollisionSample.size() + atracks.size();
      float weight = 1.;
      weight = histoReweight->GetBinContent(histoReweight->FindBin(z));

      registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z, weight);
      for (auto& track : atracks) {

        float phi = track.phis();
        o2::math_utils::bringTo02Pi(phi);

        weight = 1.; // change this !
        registry.fill(HIST("TracksEtaZvtx"), track.etas(), z, weight);
        registry.fill(HIST("TracksPhiEta"), phi, track.etas(), weight);
        registry.fill(HIST("TracksPhiZvtx"), phi, z, weight);
        registry.fill(HIST("TracksPtEta"), track.pts(), track.etas(), weight);
      }

      for (auto& track : tracks) {
        if (find(ambTrackIds.begin(), ambTrackIds.end(), track.globalIndex()) != ambTrackIds.end()) {
          continue; // this track has already been reassigned to bestcollision, don't double count
        }
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z, weight);
        registry.fill(HIST("TracksPhiEta"), phi, track.eta(), weight);
        registry.fill(HIST("TracksPtEta"), track.pt(), track.eta(), weight);
        if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
          registry.fill(HIST("TracksPhiZvtx"), phi, z, weight);
        }
      }

    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReweight, "Process reco info w reweight", false);

  using Particles = soa::Filtered<aod::McParticles>;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = (aod::mcparticle::eta < -2.8f) && (aod::mcparticle::eta > -3.2f);

  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, Particles const& particles, aod::MFTTracks const& tracks)
  {
    registry.fill(HIST("EventEfficiency"), 1.);

    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto nCharged = 0;
    for (auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nCharged++;
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());

    bool atLeastOne = false;
    int moreThanOne = 0;

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 3.);
      if (!useEvSel || (useEvSel && collision.sel8())) {
        atLeastOne = true;
        auto perCollisionSample = sample->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex());

        registry.fill(HIST("EventEfficiency"), 4.);

        registry.fill(HIST("EventsZposDiff"), collision.posZ() - mcCollision.posZ());
        if (useZDiffCut) {
          if (std::abs(collision.posZ() - mcCollision.posZ()) > maxZDiff) {
            continue;
          }
        }
        registry.fill(HIST("EventsNtrkZvtxGen"), perCollisionSample.size(), collision.posZ());
        ++moreThanOne;
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    if (moreThanOne > 1) {
      registry.fill(HIST("EventsSplitMult"), nCharged);
    }

    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }

      registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
      if (atLeastOne) {
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), mcCollision.posZ());
        registry.fill(HIST("TracksPtEtaGen"), particle.pt(), particle.eta());
      }
      registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
      registry.fill(HIST("TracksPhiZvtxGen"), particle.phi(), mcCollision.posZ());
      registry.fill(HIST("TracksPtEtaGen_t"), particle.pt(), particle.eta());
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGen, "Process generator-level info", false);

  void processGenPt(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, MFTTracksLabeled const& tracks, aod::McParticles const&)
  {
    // In the MFT the measurement of pT is not precise, so we access it by using the particle's pT instead
    if (!useEvSel || (useEvSel && collision.sel8())) {
      for (auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto particle = track.mcParticle();
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("TracksToPartPtEta"), particle.pt(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenPt, "Process particle-level info of pt", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
