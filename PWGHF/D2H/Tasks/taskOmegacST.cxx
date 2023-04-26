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

#include "TDatabasePDG.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTracking.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskOmegacST {
  using TrackedCascades = soa::Join<aod::TrackedCascades, aod::TrackedCascadeColls>;
  using TrackedV0s = soa::Join<aod::TrackedV0s, aod::TrackedV0Colls>;
  using Tracked3Bodys = soa::Join<aod::Tracked3Bodys, aod::Tracked3BodyColls>;
  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels>;

  Configurable<double> bz{"bz", 50., "magnetic field"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  OutputObj<TH1F> hDCA{"h_DCA"};
  OutputObj<TH1F> hMassOmegac{"h_MassOmegac"};
  OutputObj<TH1F> hMassOmegacMc{"h_MassOmegac_mc"};
  OutputObj<TH1F> hDCAxy{"h_DCAxy"};
  OutputObj<TH1F> hDCAz{"h_DCAz"};
  OutputObj<TH2F> hDCAVsPt{"h_DCAVsPt"};
  OutputObj<TH2F> hDCAVsR{"h_DCAVsR"};
  OutputObj<TH2F> hMassVsPt{"h_MassVsPt"};

  void init(InitContext const&)
  {
    hDCA.setObject(new TH1F("h_dca", "DCA;DCA (cm)", 200, 0., .5));
    hMassOmegac.setObject(new TH1F("h_mass_omegac", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2})", 400, 1.5, 3.));
    hMassOmegacMc.setObject(new TH1F("h_mass_omegac_mc", "inv. mass #Omega + #pi (from MC);inv. mass (GeV/#it{c}^{2})", 400, 1.5, 3.));
    // hMassOmegac.setObject(new TH1F("h_mass_omegac", "Omega_{c} mass;inv. mass (GeV/#it{c}^{2})", 200, 2.685, 2.715));
    hDCAxy.setObject(new TH1F("h_dcaxy", "DCA xy;DCA_{xy} (cm)", 200, -.5, .5));
    hDCAz.setObject(new TH1F("h_dcaz", "DCA z;DCA_{z} (cm)", 200, -.5, .5));
    hDCAVsPt.setObject(new TH2F("h_dcavspt", "DCA vs p_{T};DCA (cm);p_{T} (GeV/#it{c})", 200, -.5, .5, 200, 0., 10.));
    hDCAVsR.setObject(new TH2F("h_dcavspt", "DCA vs R;DCA (cm);R (cm)", 200, -.5, .5, 200, 0., 10.));
    hMassVsPt.setObject(new TH2F("h_massvspt", "Mass vs p_{T};Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", 200, 0., 10., 200, 0., 10.));
  }

  void process(aod::Collision const& collision,
               TrackedCascades const& trackedCascades, aod::Cascades const& cascades,
               aod::V0s const& v0s, TracksExt const& tracks, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (mRunNumber != bc.runNumber()) {
      mRunNumber = bc.runNumber();
      auto run3grp_timestamp = bc.timestamp();

      if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else {
        if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp)) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
        } else {
          LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
        }
      }
    }

    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterTrk;
    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.track_as<TracksExt>();
      auto trackParCovTrk = getTrackParCov(trackCasc);
      // trackParCovTrk.propagateToDCA(primaryVertex, bz, &impactParameterTrk);
      // gpu::gpustd::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovTrk, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &impactParameterTrk);

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();

      hDCA->Fill(TMath::Sqrt(impactParameterTrk.getR2()));
      hDCAxy->Fill(impactParameterTrk.getY());
      hDCAz->Fill(impactParameterTrk.getZ());
      hDCAVsPt->Fill(impactParameterTrk.getY(), trackCasc.pt());
      hDCAVsR->Fill(impactParameterTrk.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      hMassVsPt->Fill(trackedCascade.omegaMass(), trackCasc.pt());

      if (!ptrack.has_mcParticle() || !ntrack.has_mcParticle() || !bachelor.has_mcParticle())
        continue;

      LOGF(debug, "ptrack (id: %d, pdg: %d) has mother %d", ptrack.mcParticleId(),
           ptrack.mcParticle().pdgCode(), ptrack.mcParticle().has_mothers() ? ptrack.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "ntrack (id: %d, pdg: %d) has mother %d", ntrack.mcParticleId(),
           ntrack.mcParticle().pdgCode(), ntrack.mcParticle().has_mothers() ? ntrack.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          const auto mother = v0part.mothers_as<aod::McParticles>()[0];
          const auto pdgCode = mother.pdgCode();
          LOG(debug) << "cascade with PDG code: " << pdgCode;
          if (std::abs(pdgCode) == 3334) {
            LOG(debug) << "found Omega, looking for pions";
            std::array<double, 2> masses{RecoDecay::getMassPDG(3334), RecoDecay::getMassPDG(211)};
            std::array<std::array<float, 3>, 2> momenta;
            trackParCovTrk.getPxPyPzGlo(momenta[0]);
            // momenta[0] = {trackCasc.px(), trackCasc.py(), trackCasc.pz()};
            for (const auto& track : tracks) {
              if (!track.has_mcParticle())
                continue;
              const auto mcpart = track.mcParticle();
              if (mcpart.pdgCode() == std::copysign(211, pdgCode)) {
                LOGF(debug, "combining Omega with pion %d", track.globalIndex());
                auto trackParCovPion = getTrackParCov(track);
                o2::dataformats::DCA impactParameterPion;
                // trackParCovPion.propagateToDCA(primaryVertex, bz, &impactParameterPion);
                o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPion, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &impactParameterPion);
                trackParCovPion.getPxPyPzGlo(momenta[1]);
                // momenta[1] = {track.px(), track.py(), track.pz()};
                hMassOmegac->Fill(RecoDecay::m(momenta, masses));

                // MC-based mass
                momenta[0] = {mother.px(), mother.py(), mother.pz()};
                momenta[1] = {mcpart.px(), mcpart.py(), mcpart.pz()};
                hMassOmegacMc->Fill(RecoDecay::m2(momenta, masses));
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskOmegacST>(cfgc)};
}
