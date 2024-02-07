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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct StrangenessTrackingQATask {
  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels>;

  Configurable<double> bz{"bz", 50., "magnetic field"};

  OutputObj<TH1F> hDCA{"h_DCA"};
  OutputObj<TH1F> hDCAxy{"h_DCAxy"};
  OutputObj<TH1F> hDCAz{"h_DCAz"};
  OutputObj<TH2F> hDCAVsPt{"h_DCAVsPt"};
  OutputObj<TH2F> hDCAVsR{"h_DCAVsR"};
  OutputObj<TH2F> hMassVsPt{"h_MassVsPt"};
  OutputObj<TH2F> hBuilderMassVsPt{"h_BuilderMassVsPt"};
  OutputObj<TH2F> hMassVsMass{"h_MassVsMass"};

  template <typename T>
  float QuadraticSum(T x, T y)
  {
    return TMath::Sqrt(x * x + y * y);
  }

  void init(InitContext const&)
  {
    hDCA.setObject(new TH1F("h_dca", "DCA;DCA (cm)", 200, 0., .5));
    hDCAxy.setObject(new TH1F("h_dcaxy", "DCA xy;DCA_{xy} (cm)", 200, -.5, .5));
    hDCAz.setObject(new TH1F("h_dcaz", "DCA z;DCA_{z} (cm)", 200, -.5, .5));
    hDCAVsPt.setObject(new TH2F("h_dcavspt", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", 200, -.5, .5, 200, 0., 10.));
    hDCAVsR.setObject(new TH2F("h_dcavspt", "DCA vs R;DCA (cm);R (cm)", 200, -.5, .5, 200, 0., 10.));
    hMassVsPt.setObject(new TH2F("h_massvspt", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", 200, 0., 10., 200, 0., 10.));
    hBuilderMassVsPt.setObject(new TH2F("h_buildermassvspt", "Mass (from builder) vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", 200, 0., 10., 200, 0., 10.));
    hMassVsMass.setObject(new TH2F("h_massvsmass", "Mass vs mass;Mass (GeV/#it{c}^{2});Mass (GeV/#it{c}^{2})", 200, 0., 10., 200, 0., 10.));
  }

  void processTrackedCascades(aod::Collision const& collision,
                              aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                              aod::V0s const& v0s, TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& trackedCascade : trackedCascades) {
      const auto track = trackedCascade.track_as<TracksExt>();
      auto trackCovTrk = getTrackParCov(track);
      auto primaryVertex = getPrimaryVertex(collision);
      // auto covMatrixPV = primaryVertex.getCov();
      o2::dataformats::DCA impactParameterTrk;
      trackCovTrk.propagateToDCA(primaryVertex, bz, &impactParameterTrk);

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();

      hDCA->Fill(TMath::Sqrt(impactParameterTrk.getR2()));
      hDCAxy->Fill(impactParameterTrk.getY());
      hDCAz->Fill(impactParameterTrk.getZ());
      hDCAVsPt->Fill(impactParameterTrk.getY(), track.pt());
      hDCAVsR->Fill(impactParameterTrk.getY(), QuadraticSum(track.x(), track.y()));
      hMassVsPt->Fill(trackedCascade.omegaMass(), track.pt());

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
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessTrackingQATask, processTrackedCascades, "process cascades from strangeness tracking", true);

  void processCascades(aod::Collision const& collision, aod::TrackedCascades const& trackedCascades, aod::Cascades const& cascades, aod::V0s const& v0s,
                       soa::Join<aod::TraCascDatas, aod::McTraCascLabels> const& trackedcascdata, TracksExt const& tracks, aod::McParticles const& mcParticles)
  {
    for (const auto& trackedCascadeData : trackedcascdata) {
      hBuilderMassVsPt->Fill(trackedCascadeData.mOmega(), trackedCascadeData.pt());
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
              hMassVsMass->Fill(trackedCascade.omegaMass(), trackedCascadeData.mOmega());
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
