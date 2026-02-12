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
//
// Analysis task to produce resolution mapfor electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

// #include "TGeoGlobalMagField.h"

#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

struct testBremsstrahlung {
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float d_bz = 0;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  ~testBremsstrahlung() {}

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    mRunNumber = 0;
    d_bz = 0;

    if (doprocessGen) {
      registry.add("Event/hGenID", "generator ID;generator ID;Number of mc collisions", kTH1F, {{7, -1.5, 5.5}}, false);
      registry.add("Photon/hXY", "prod. vtx. of #gamma^{bremss};X (cm);Y (cm)", kTH2D, {{400, -100, +100}, {400, -100, +100}}, false);
      registry.add("Photon/hs", "kinetic var. at prod. vtx.;p_{T,#gamma}^{bremss} (GeV/c);#eta_{#gamma}^{bremss};#varphi_{#gamma}^{bremss} (rad.);", kTHnSparseD, {{100, 0, 1}, {40, -1, +1}, {72, 0, 2 * M_PI}}, false);
      registry.add("Photon/hsDelta", "hsDelta;(p_{T,e}^{after} - p_{T,e}^{before})/p_{T,e}^{before};p_{T,e}^{before} (GeV/c);r_{xy}^{bremss} (cm);#varphi^{bremss} (rad.);z^{bremss} (cm);", kTHnSparseD, {{100, -1, 0}, {100, 0, 1}, {200, 0, 100}, {72, 0, 2 * M_PI}, {200, -100, +100}}, false);
    }
  }

  // void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  // {
  //   if (mRunNumber == bc.runNumber()) {
  //     return;
  //   }

  //   // load matLUT for this timestamp
  //   if (!lut) {
  //     LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
  //     lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
  //   } else {
  //     LOG(info) << "Material look-up table already in place. Not reloading.";
  //   }

  //   // In case override, don't proceed, please - no CCDB access required
  //   if (d_bz_input > -990) {
  //     d_bz = d_bz_input;
  //     o2::parameters::GRPMagField grpmag;
  //     if (std::fabs(d_bz) > 1e-5) {
  //       grpmag.setL3Current(30000.f / (d_bz / 5.0f));
  //     }
  //     o2::base::Propagator::initFieldFromGRP(&grpmag);
  //     o2::base::Propagator::Instance()->setMatLUT(lut);
  //     mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
  //     mRunNumber = bc.runNumber();

  //     if (!o2::base::GeometryManager::isGeometryLoaded()) {
  //       ccdb->get<TGeoManager>(geoPath);
  //     }
  //   }

  //   auto run3grp_timestamp = bc.timestamp();
  //   o2::parameters::GRPObject* grpo = 0x0;
  //   o2::parameters::GRPMagField* grpmag = 0x0;
  //   if (!skipGRPOquery) {
  //     grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
  //   }
  //   if (grpo) {
  //     o2::base::Propagator::initFieldFromGRP(grpo);
  //     o2::base::Propagator::Instance()->setMatLUT(lut);
  //     mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
  //     // Fetch magnetic field from ccdb for current collision
  //     d_bz = grpo->getNominalL3Field();
  //     LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
  //   } else {
  //     grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
  //     if (!grpmag) {
  //       LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
  //     }
  //     o2::base::Propagator::initFieldFromGRP(grpmag);
  //     o2::base::Propagator::Instance()->setMatLUT(lut);
  //     mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());

  //     // Fetch magnetic field from ccdb for current collision
  //     d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
  //     LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
  //   }

  //   mRunNumber = bc.runNumber();
  // }

  SliceCache cache;
  Preslice<aod::TracksIU> perCollision_mid = o2::aod::track::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  using MyCollisions = Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollision = MyCollisions::iterator;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
  using MyTrack = MyTracks::iterator;

  void processGen(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcCollisions) {
      registry.fill(HIST("Event/hGenID"), mcCollision.getSubGeneratorId());

      auto mcParticles_per_collision = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
      // auto photons_per_collision = photons->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

      for (const auto& mcParticle : mcParticles_per_collision) {
        if (std::abs(mcParticle.pdgCode()) != 22) {
          continue;
        }
        if (!mcParticle.has_mothers()) { // select seconday mcParticles
          continue;
        }
        if (mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator()) { // select seconday mcParticles
          continue;
        }

        auto mother = mcParticle.template mothers_as<aod::McParticles>()[0]; // to be primary electron
        if (mother.eta() < -0.9 || 0.9 < mother.eta()) {
          continue;
        }
        if (std::abs(mother.pdgCode()) != 11) { // select bremsstrahlung photon from primary electron
          continue;
        }
        if (!(mother.isPhysicalPrimary() || mother.producedByGenerator())) { // select bremsstrahlung photon from primary electron
          continue;
        }

        registry.fill(HIST("Photon/hXY"), mcParticle.vx(), mcParticle.vy());
        registry.fill(HIST("Photon/hs"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi() > 0 ? mcParticle.phi() : mcParticle.phi() + 2 * M_PI);
        // int ndau = mother.daughtersIds()[1] - mother.daughtersIds()[0] + 1;
        // LOGF(info, "ndau = %d", ndau);

        for (const auto& daughter : mother.daughters_as<aod::McParticles>()) {
          // LOGF(info, "mother.globalIndex() = %d, daughter.globalIndex() = %d, daughter.pdgCode() = %d", mother.globalIndex(), daughter.globalIndex(), daughter.pdgCode());
          if (std::abs(daughter.pdgCode()) == 11 && (!daughter.isPhysicalPrimary() && !daughter.producedByGenerator())) {
            float reldpt = (daughter.pt() - mother.pt()) / mother.pt();
            float rxy = std::sqrt(std::pow(daughter.vx(), 2) + std::pow(daughter.vy(), 2));
            float phi_bremss = RecoDecay::phi(daughter.vy(), daughter.vx());
            o2::math_utils::bringTo02Pi(phi_bremss);
            registry.fill(HIST("Photon/hsDelta"), reldpt, mother.pt(), rxy, phi_bremss, daughter.vz());
          }
        }

      } // end of mc particle per mc collision

    } // end of mc collision loop
  }
  PROCESS_SWITCH(testBremsstrahlung, processGen, "process generated info", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<testBremsstrahlung>(cfgc, TaskName{"test-bremsstrahlung"})};
}
