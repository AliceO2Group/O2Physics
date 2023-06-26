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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"



using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct StrangenessTrackingQATask {
  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels>;

  Configurable<double> bz{"bz", -50., "magnetic field"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0.f;


  HistogramRegistry registry{
    "registry",
    {
      {"h_dca", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"h_dcaxy", "DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_dcaz", "DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"h_bachdcaxy", "Bachelor DCA xy;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_bachdcaz", "Bachelor DCA z;DCA_{z} (cm)", {HistType::kTH1D, {{200, -1., 1.}}}},
      {"h_dcavspt", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavspt", "Bachelor DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_bachdcavsr", "Bachelor DCA vs R (cm);DCA (cm);R (cm)", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 30.}}}},
      {"h_ntrackdcavspt", "N track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_ptrackdcavspt", "P track DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, -1., 1.}, {200, 0., 10.}}}},
      {"h_dcavsr", "DCA vs R;DCA (cm);R (cm)", {HistType::kTH2D, {{200, -.5, .5}, {200, 0., 10.}}}},
      {"h_massvspt", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_buildermassvspt", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D,  {{200, 0., 10.}, {200, 0., 10.}}}},
      {"h_massvsmass", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
    }};


  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, run3grp_timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, run3grp_timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    } else {
      LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPpath << " of object GRPMagField and " << cfgGRPmagPath << " of object GRPObject for timestamp " << run3grp_timestamp;
    }
  }



  void init(InitContext const&)
  {
    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }
  }

  void processTrackedCascades(aod::Collision const& collision,
                              aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                              aod::V0s const& v0s, TracksExt const& tracks, aod::McParticles const& mcParticles,
                              aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);

    // gpu::gpustd::array<float, 2> dca{-999.f, -999.f};
    const auto primaryVertex = getPrimaryVertex(collision);
    for (const auto& trackedCascade : trackedCascades) {
      const auto track = trackedCascade.track_as<TracksExt>();
      auto trackCovTrk = getTrackParCov(track);
      o2::dataformats::DCA impactParameterTrk;
      // trackCovTrk.propagateToDCA(primaryVertex, bz, &impactParameterTrk);
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovTrk, mBz, 2.f, matCorr, &impactParameterTrk)){
        registry.fill(HIST("h_dca"), TMath::Sqrt(impactParameterTrk.getR2()));
        registry.fill(HIST("h_dcaxy"), impactParameterTrk.getY());
        registry.fill(HIST("h_dcaz"), impactParameterTrk.getZ());
        registry.fill(HIST("h_dcavspt"), impactParameterTrk.getY(), track.pt());
        registry.fill(HIST("h_dcavsr"), impactParameterTrk.getY(), std::hypot(track.x(),track.y() ));
        registry.fill(HIST("h_massvspt"), trackedCascade.omegaMass(), track.pt());
      }

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();


      LOGF(info, "ptrack (id: %d, pdg: %d) has mother %d", ptrack.mcParticleId(),
           ptrack.mcParticle().pdgCode(), ptrack.mcParticle().has_mothers() ? ptrack.mcParticle().mothersIds()[0] : -1);
      LOGF(info, "ntrack (id: %d, pdg: %d) has mother %d", ntrack.mcParticleId(),
           ntrack.mcParticle().pdgCode(), ntrack.mcParticle().has_mothers() ? ntrack.mcParticle().mothersIds()[0] : -1);

      LOG(info) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        LOG(info) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          LOG(info) << "cascade with PDG code: " << v0part.mothers_as<aod::McParticles>()[0].pdgCode();
        } else {
          continue;
        }
      }
      auto trackCovBach = getTrackParCov(bachelor);
      o2::dataformats::DCA impactParameterBach;
      // trackCovBach.propagateToDCA(primaryVertex, bz, &impactParameterBach);
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovBach, mBz, 2.f, matCorr, &impactParameterBach)){
        registry.fill(HIST("h_bachdcaxy"), impactParameterBach.getY());
        registry.fill(HIST("h_bachdcaz"), impactParameterBach.getZ());
        registry.fill(HIST("h_bachdcavspt"), impactParameterBach.getY(), bachelor.pt());
        registry.fill(HIST("h_bachdcavsr"), impactParameterBach.getY(),std::hypot(trackedCascade.decayX(),trackedCascade.decayY()));
      }


      auto trackCovNtrack = getTrackParCov(ntrack);
      o2::dataformats::DCA impactParameterNtrack;
      // trackCovNtrack.propagateToDCA(primaryVertex, bz, &impactParameterNtrack);
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovNtrack, mBz, 2.f, matCorr, &impactParameterNtrack)){
        registry.fill(HIST("h_ntrackdcavspt"),impactParameterNtrack.getY(), ntrack.pt());
      }
      

      auto trackCovPtrack = getTrackParCov(ptrack);
      o2::dataformats::DCA impactParameterPtrack;
      // trackCovPtrack.propagateToDCA(primaryVertex, bz, &impactParameterPtrack);
      if (o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackCovPtrack, mBz, 2.f, matCorr, &impactParameterPtrack)){
        registry.fill(HIST("h_ptrackdcavspt"),impactParameterPtrack.getY(), ptrack.pt());
      }

    }
  }
  PROCESS_SWITCH(StrangenessTrackingQATask, processTrackedCascades, "process cascades from strangeness tracking", true);

  void processCascades(aod::Collision const& collision, aod::TrackedCascades const& trackedCascades, aod::Cascades const& cascades, aod::V0s const& v0s,
                       soa::Join<aod::TraCascDatas, aod::McTraCascLabels> const& trackedcascdata, TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& trackedCascadeData : trackedcascdata) {
      registry.fill(HIST("h_buildermassvspt"),trackedCascadeData.mOmega(), trackedCascadeData.pt());
    }

    for (const auto& trackedCascade : trackedCascades) {
      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          int mcid = v0part.mothersIds()[0];
          for (const auto& trackedCascadeData : trackedcascdata) {
            if (trackedCascadeData.mcParticleId() == mcid) {
              registry.fill(HIST("h_massvsmass"),trackedCascade.omegaMass(), trackedCascadeData.mOmega());

              break;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessTrackingQATask, processCascades, "process cascades from builder", true);

  void processTrackedV0s(aod::Collision const& collision,
                         aod::AssignedTrackedV0s const& trackedV0s, aod::V0s const& v0s,
                         TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& trackedV0 : trackedV0s) {
      const auto& v0 = trackedV0.v0();
      v0.posTrack();
      v0.negTrack();
    }
  }
  PROCESS_SWITCH(StrangenessTrackingQATask, processTrackedV0s, "process tracked V0s", true);

  void processTracked3Bodys(aod::Collision const& collision,
                            aod::AssignedTracked3Bodys const& tracked3Bodys, aod::Decay3Bodys const& decay3Bodys,
                            TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& tracked3Body : tracked3Bodys) {
      tracked3Body.itsTrack();
      const auto& decay3Body = tracked3Body.decay3Body();
      decay3Body.track0();
      decay3Body.track1();
      decay3Body.track2();
    }
  }
  PROCESS_SWITCH(StrangenessTrackingQATask, processTracked3Bodys, "process tracked 3 body decays", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessTrackingQATask>(cfgc)};
}
