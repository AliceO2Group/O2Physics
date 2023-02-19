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
// ========================
//
// This code will create data table for inputs to machine learning to classify V0s.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include <map>
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using FullTrackExt = FullTracksExt::iterator;

using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;
using FullTrackExtMC = FullTracksExtMC::iterator;

namespace o2::aod
{

namespace mycollision // reconstructed collision information
{
}
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, bc::RunNumber, collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib, evsel::Sel8,
                  mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV);
using MyCollision = MyCollisions::iterator;

namespace mymccollision // generated collision information
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, Mycollision); //!
}
DECLARE_SOA_TABLE(MyMCCollisions, "AOD", "MYMCCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, mymccollision::MyCollisionId,
                  mccollision::GeneratorsID, mccollision::PosX, mccollision::PosY, mccollision::PosZ);
using MyCollision = MyCollisions::iterator;

namespace mytrack
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
DECLARE_SOA_COLUMN(Sign, sign, int);            //!
} // namespace mytrack

// reconstructed track information
DECLARE_SOA_TABLE(MyTracks, "AOD", "MYTRACK", //!
                  o2::soa::Index<>, mytrack::CollisionId, mytrack::Sign,
                  track::Pt, track::Eta, track::Phi,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap,
                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);

// iterators
using MyTrack = MyTracks::iterator;

namespace mymctrack
{
DECLARE_SOA_INDEX_COLUMN(MyTrack, mytrack);                     //!
DECLARE_SOA_COLUMN(MotherPdgCode, motherpdgCode, int);          //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //!
} // namespace mymctrack

// generated track information
DECLARE_SOA_TABLE(MyMCTracks, "AOD", "MYMCTRACK", //!
                  o2::soa::Index<>, mymctrack::MyTrackId,
                  mcparticle::Pt, mcparticle::Eta, mcparticle::Phi,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,
                  mcparticle::PdgCode, mymctrack::IsPhysicalPrimary, mymctrack::MotherPdgCode);
using MyMCTrack = MyMCTracks::iterator;

} // namespace o2::aod

struct TreeCreatorElectronML {
  Produces<o2::aod::MyCollisions> mycollision;
  Produces<o2::aod::MyMCCollisions> mymccollision;
  Produces<o2::aod::MyTracks> mytrack;
  Produces<o2::aod::MyMCTracks> mymctrack;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
    },
  };

  // Configurables
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};

  void init(InitContext& context) {}

  Filter trackFilter = nabs(o2::aod::track::eta) < maxeta && o2::aod::track::tpcChi2NCl < maxchi2tpc;
  using MyFilteredTracksMC = soa::Filtered<FullTracksExtMC>;
  Preslice<MyFilteredTracksMC> perCollision = aod::track::collisionId;

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const& particlesMC, aod::McCollisions const&)
  {
    for (auto& collision : collisions) {
      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }

      registry.fill(HIST("hEventCounter"), 1.0); // all

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      mycollision(bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.sel8(),
                  collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                  collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV());

      auto mccollision = collision.mcCollision();
      mymccollision(mycollision.lastIndex(), mccollision.generatorsID(), mccollision.posX(), mccollision.posY(), mccollision.posZ());

      auto tracks_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto& track : tracks_coll) {

        if (track.tpcNClsCrossedRows() < mincrossedrows) {
          continue;
        }

        mytrack(mycollision.lastIndex(),
                track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(),
                track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap());

        auto mcparticle = track.mcParticle_as<aod::McParticles>();

        if (mcparticle.has_mothers()) {
          auto fm = mcparticle.mothers_as<aod::McParticles>().front(); // first mother
          mymctrack(mytrack.lastIndex(),
                    mcparticle.pt(), mcparticle.eta(), mcparticle.phi(),
                    mcparticle.vx(), mcparticle.vy(), mcparticle.vz(),
                    mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), fm.pdgCode());
        } else {
          mymctrack(mytrack.lastIndex(),
                    mcparticle.pt(), mcparticle.eta(), mcparticle.phi(),
                    mcparticle.vx(), mcparticle.vy(), mcparticle.vz(),
                    mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), 0);
        }

      } // end of track loop
    }   // end of collision loop
  }     // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronML>(cfgc, TaskName{"tree-creator-ele-ml"})};
}
