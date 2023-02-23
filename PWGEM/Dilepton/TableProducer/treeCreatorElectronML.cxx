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
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using FullTrackExt = FullTracksExt::iterator;

using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;
using FullTrackExtMC = FullTracksExtMC::iterator;

namespace o2::aod
{

namespace mycollision // reconstructed collision information
{
DECLARE_SOA_COLUMN(MCPosX, mcposX, float); //!
DECLARE_SOA_COLUMN(MCPosY, mcposY, float); //!
DECLARE_SOA_COLUMN(MCPosZ, mcposZ, float); //!
} // namespace mycollision
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, bc::GlobalBC, bc::RunNumber, collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib, evsel::Sel8,
                  mccollision::GeneratorsID, mycollision::MCPosX, mycollision::MCPosY, mycollision::MCPosZ,
                  mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV, mult::MultNTracksPVeta1);
using MyCollision = MyCollisions::iterator;

namespace mytrack
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);              //!
DECLARE_SOA_COLUMN(Sign, sign, int);                             //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);             //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int); //!
DECLARE_SOA_COLUMN(MCPt, mcpt, float);                           //!
DECLARE_SOA_COLUMN(MCEta, mceta, float);                         //!
DECLARE_SOA_COLUMN(MCPhi, mcphi, float);                         //!
DECLARE_SOA_COLUMN(MCVx, mcvx, float);                           //!
DECLARE_SOA_COLUMN(MCVy, mcvy, float);                           //!
DECLARE_SOA_COLUMN(MCVz, mcvz, float);                           //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);  //!
DECLARE_SOA_COLUMN(MotherPdgCode, motherpdgCode, int);           //!
DECLARE_SOA_COLUMN(GrandMotherPdgCode, grandmotherpdgCode, int); //!
} // namespace mytrack

// reconstructed track information
DECLARE_SOA_TABLE(MyTracks, "AOD", "MYTRACK", //!
                  o2::soa::Index<>, mytrack::MyCollisionId, mytrack::Sign,
                  track::Pt, track::Eta, track::Phi,
                  track::DcaXY, track::DcaZ,
                  // track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsFindable, mytrack::TPCNClsFound, mytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap,
                  mytrack::MCPt, mytrack::MCEta, mytrack::MCPhi,
                  mytrack::MCVx, mytrack::MCVy, mytrack::MCVz,
                  mcparticle::PdgCode, mytrack::IsPhysicalPrimary, mytrack::MotherPdgCode, mytrack::GrandMotherPdgCode,
                  // dynamic column
                  // track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);

// iterators
using MyTrack = MyTracks::iterator;

} // namespace o2::aod

struct TreeCreatorElectronML {
  Produces<o2::aod::MyCollisions> mycollision;
  Produces<o2::aod::MyTracks> mytrack;

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
      auto mccollision = collision.mcCollision();
      mycollision(bc.globalBC(), bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.sel8(),
                  mccollision.generatorsID(), mccollision.posX(), mccollision.posY(), mccollision.posZ(),
                  collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                  collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1());

      auto tracks_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto& track : tracks_coll) {

        if (track.tpcNClsCrossedRows() < mincrossedrows) {
          continue;
        }
        int mother_pdg = 0;
        int grand_mother_pdg = 0;
        auto mcparticle = track.mcParticle_as<aod::McParticles>();

        if (mcparticle.has_mothers()) {
          auto mother = mcparticle.mothers_as<aod::McParticles>().front(); // first mother
          mother_pdg = mother.pdgCode();
          if (mother.has_mothers()) {
            auto grand_mother = mother.mothers_as<aod::McParticles>().front(); // first mother of mother
            grand_mother_pdg = grand_mother.pdgCode();
          }
        }

        mytrack(mycollision.lastIndex(),
                track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(),
                track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(),
                mcparticle.pt(), mcparticle.eta(), mcparticle.phi(),
                mcparticle.vx(), mcparticle.vy(), mcparticle.vz(),
                mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), mother_pdg, grand_mother_pdg);

      } // end of track loop
    }   // end of collision loop
  }     // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronML>(cfgc, TaskName{"tree-creator-ele-ml"})};
}
