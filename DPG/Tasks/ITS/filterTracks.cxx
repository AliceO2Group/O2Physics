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
/// \author Andrea Rossi <andrea.rossi@cern.ch>

/// \brief Simple task to filter tracks and save infos to trees for DCA-related studies (alignment, HF-related issues)

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;
using namespace o2::aod::mctracklabel;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace filterTracks
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision);              //! Collision index
DECLARE_SOA_COLUMN(IsInsideBeamPipe, isInsideBeamPipe, int); //! is within beam pipe
DECLARE_SOA_COLUMN(Pt, pt, float);                           //! track pt
DECLARE_SOA_COLUMN(Px, px, float);                           //! track px
DECLARE_SOA_COLUMN(Py, py, float);                           //! track py
DECLARE_SOA_COLUMN(Pz, pz, float);                           //! track pz
// DECLARE_SOA_COLUMN(Eta, eta, float);                                //! track eta
// DECLARE_SOA_COLUMN(X, x, float);                          //! track x position at the DCA to the primary vertex
// DECLARE_SOA_COLUMN(Y, y, float);                          //! track y position at the DCA to the primary vertex
// DECLARE_SOA_COLUMN(Z, z, float);                          //! track z position at the DCA to the primary vertex
// DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                          //! track distance of closest approach at the primary vertex: in xy plane
// DECLARE_SOA_COLUMN(DcaZ, dcaz, float);               //! track distance of closest approach at the primary vertex: along z (beam line) direction
DECLARE_SOA_COLUMN(Charge, charge, int);             //! track sign, not really charge
DECLARE_SOA_COLUMN(NsigmaTPCpi, nsigmaTPCpi, float); //! TPC nsigma w.r.t. pion mass hypothesis
DECLARE_SOA_COLUMN(NsigmaTPCka, nsigmaTPCka, float); //! TPC nsigma w.r.t. kaon mass hypothesis
DECLARE_SOA_COLUMN(NsigmaTPCpr, nsigmaTPCpr, float); //! TPC nsigma w.r.t. proton mass hypothesis
DECLARE_SOA_COLUMN(NsigmaTOFpi, nsigmaTOFpi, float); //! TOF nsigma w.r.t. pion mass hypothesis
DECLARE_SOA_COLUMN(NsigmaTOFka, nsigmaTOFka, float); //! TOF nsigma w.r.t. kaon mass hypothesis
DECLARE_SOA_COLUMN(NsigmaTOFpr, nsigmaTOFpr, float); //! TOF nsigma w.r.t. proton mass hypothesis
DECLARE_SOA_COLUMN(TpcNCluster, tpcNCluster, int);   //! TOF nsigma w.r.t. proton mass hypothesis

///// MC INFO
DECLARE_SOA_COLUMN(MainHfMotherPdgCode, mainMotherPdgCode, int);                   //! mother pdg code for particles coming from HF, skipping intermediate resonance states. Not trustable when mother is not HF. Not suited for Sc->Lc decays, since Sc are never pointed to
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);                    //! is phyiscal primary according to ALICE definition
DECLARE_SOA_COLUMN(MainBeautyAncestorPdgCode, mainBeautyAncestorPdgCode, int);     //! pdgcode of beauty particle when this is an ancestor, otherwise -1
DECLARE_SOA_COLUMN(MainMotherOrigIndex, mainMotherOrigIndex, int);                 //! original index in MCParticle tree of main mother: needed when checking if particles come from same mother
DECLARE_SOA_COLUMN(MainMotherNfinalStateDaught, mainMotherNfinalStateDaught, int); //! number of final state (we consider only pions, kaons, muons, electrons, protons) daughter in main mother decay. To be noted that this is computed only for decays of particles of interest (D0, Lc, K0s). If the sign is negative, it means that the decay is not in one of the desired channels (K0s->pi pi, Lc->pKpi, D0->K-pi+)

DECLARE_SOA_COLUMN(MainMotherPt, mainMotherPt, float);                 //! original index in MCParticle tree of main mother: needed when chekcing if particles come from same mother
DECLARE_SOA_COLUMN(MainMotherY, mainMotherY, float);                   //! original index in MCParticle tree of main mother: needed when chekcing if particles come from same mother
DECLARE_SOA_COLUMN(MainBeautyAncestorPt, mainBeautyAncestorPt, float); //! original index in MCParticle tree of main mother: needed when chekcing if particles come from same mother
DECLARE_SOA_COLUMN(MainBeautyAncestorY, mainBeautyAncestorY, float);   //! original index in MCParticle tree of main mother: needed when chekcing if particles come from same mother
DECLARE_SOA_COLUMN(MaxEtaDaughter, maxEtaDaughter, float);             //! max (abs) eta of daughter particles, needed to reproduce acceptance cut
} // namespace filterTracks
DECLARE_SOA_TABLE(FilterTrack, "AOD", "FILTERTRACK",
                  o2::aod::track::CollisionId,
                  aod::filterTracks::IsInsideBeamPipe,
                  o2::aod::track::TrackType,
                  o2::aod::track::X,
                  o2::aod::track::Alpha,
                  o2::aod::track::Y,
                  o2::aod::track::Z,
                  o2::aod::track::Snp,
                  o2::aod::track::Tgl,
                  o2::aod::track::Signed1Pt);
DECLARE_SOA_TABLE(FilterTrackExtr, "AOD", "FILTERTRACKEXTR",
                  // aod::filterTracks::Px,aod::filterTracks::Py, aod::filterTracks::Pz,
                  aod::filterTracks::Pt, o2::aod::track::Eta,
                  o2::aod::filterTracks::Charge,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  o2::aod::track::SigmaDcaXY2,
                  o2::aod::track::SigmaDcaZ2,
                  aod::filterTracks::NsigmaTPCpi, aod::filterTracks::NsigmaTPCka, aod::filterTracks::NsigmaTPCpr,
                  aod::filterTracks::NsigmaTOFpi, aod::filterTracks::NsigmaTOFka, aod::filterTracks::NsigmaTOFpr);
DECLARE_SOA_TABLE(FiltTracExtDet, "AOD", "FILTTRACEXTDET",
                  o2::aod::track::ITSClusterSizes,
                  o2::aod::track::ITSChi2NCl,
                  o2::aod::track::TPCChi2NCl,
                  aod::filterTracks::TpcNCluster,
                  o2::aod::track::TrackTime);
DECLARE_SOA_TABLE(FilterTrackMC, "AOD", "FILTERTRACKMC",
                  // aod::filterTracks::Px,aod::filterTracks::Py, aod::filterTracks::Pz,
                  o2::aod::mcparticle::PdgCode,
                  o2::aod::filterTracks::IsPhysicalPrimary,
                  o2::aod::filterTracks::MainHfMotherPdgCode,
                  o2::aod::filterTracks::MainBeautyAncestorPdgCode,
                  o2::aod::filterTracks::MainMotherOrigIndex,
                  o2::aod::filterTracks::MainMotherNfinalStateDaught,
                  o2::aod::filterTracks::MainMotherPt,
                  o2::aod::filterTracks::MainMotherY,
                  o2::aod::filterTracks::MainBeautyAncestorPt,
                  o2::aod::filterTracks::MainBeautyAncestorY);
DECLARE_SOA_TABLE(GenParticles, "AOD", "GENPARTICLES",
                  // aod::filterTracks::Px,aod::filterTracks::Py, aod::filterTracks::Pz,
                  o2::aod::mcparticle::PdgCode,
                  o2::aod::mcparticle::McCollisionId,
                  o2::aod::filterTracks::MainBeautyAncestorPdgCode,
                  o2::aod::filterTracks::MainMotherPt,
                  o2::aod::filterTracks::MainMotherY,
                  o2::aod::filterTracks::MaxEtaDaughter,
                  o2::aod::filterTracks::MainBeautyAncestorPt,
                  o2::aod::filterTracks::MainBeautyAncestorY);
} // namespace o2::aod

struct FilterTracks {

  Produces<aod::FilterTrackExtr> filteredTracksTableExtra;
  Produces<aod::FilterTrack> filteredTracksTable;
  Produces<aod::FiltTracExtDet> filteredTracksTableExtraDet;
  Produces<aod::FilterTrackMC> filteredTracksMC;
  Produces<aod::GenParticles> selectedGenParticles;

  //  Configurable<int> dummy{"dummy", 0, "dummy"};
  Configurable<float> minTrackPt{"minTrackPt", 0.25, "min track pt"};
  Configurable<float> trackDcaXyMax{"trackDcaXyMax", 0.5, "max track pt"};
  Configurable<int> trackPtSampling{"trackPtSampling", 0, "track sampling mode"};
  Configurable<float> trackPtWeightLowPt{"trackPtWeightLowPt", 0.01f, "trackPtWeightLowPt"};
  Configurable<float> trackPtWeightMidPt{"trackPtWeightMidPt", 0.10f, "trackPtWeightMidPt"};

  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && aod::track::pt > minTrackPt&& nabs(aod::track::dcaXY) < trackDcaXyMax;
  using TracksWithSelAndDca = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using TracksWithSelAndDcaMc = soa::Join<TracksWithSelAndDca, aod::McTrackLabels>;
  Partition<soa::Filtered<TracksWithSelAndDca>> lowPtTracks = aod::track::pt < 2.f && (nabs(aod::track::pt * 10000.f - nround(aod::track::pt * 10000.f)) < trackPtWeightLowPt * 2.f);
  Partition<soa::Filtered<TracksWithSelAndDca>> midPtTracks = aod::track::pt > 2.f && aod::track::pt < 5.f && (nabs(aod::track::pt * 10000.f - nround(aod::track::pt * 10000.f)) < trackPtWeightMidPt * 2.f);
  Partition<soa::Filtered<TracksWithSelAndDca>> highPtTracks = aod::track::pt > 5.f;

  Partition<soa::Filtered<TracksWithSelAndDcaMc>> lowPtTracksMC = aod::track::pt < 2.f && (nabs(aod::track::pt * 10000.f - nround(aod::track::pt * 10000.f)) < trackPtWeightLowPt * 2.f);
  Partition<soa::Filtered<TracksWithSelAndDcaMc>> midPtTracksMC = aod::track::pt > 2.f && aod::track::pt < 5.f && (nabs(aod::track::pt * 10000.f - nround(aod::track::pt * 10000.f)) < trackPtWeightMidPt * 2.f);
  Partition<soa::Filtered<TracksWithSelAndDcaMc>> highPtTracksMC = aod::track::pt > 5.f;

  std::array<int, 3> pdgSignalParticleArray = {310, 421, 4122}; // K0s, D0 and Lc
  std::array<int, 3> pdgDecayLc = {2212, -321, 211};
  std::array<int, 2> pdgDecayDzero = {-321, 211};
  std::array<int, 2> pdgDecayKzero = {-211, 211};

  void init(InitContext&)
  {
  }

  void fillTableData(auto track)
  {

    filteredTracksTableExtra(track.pt(), track.eta(), track.sign(), track.dcaXY(), track.dcaZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
    filteredTracksTable(track.collisionId(), track.isWithinBeamPipe() ? 1 : 0, track.trackType(), track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt());

    filteredTracksTableExtraDet(track.itsClusterSizes(), track.itsChi2NCl(), track.tpcChi2NCl(), track.tpcNClsFound(), track.trackTime());
  }

  void fillTableDataMC(auto track, aod::McParticles const& mcParticles)
  {

    fillTableData(track);
    bool has_MCparticle = track.has_mcParticle();
    if (has_MCparticle) {
      /// the track is not fake

      // check whether the particle comes from a charm or beauty hadron and store its index

      auto mcparticle = track.mcParticle();
      int pdgParticleMother = 0;
      for (int iSignPart = 0; iSignPart < 3; iSignPart++) {
        pdgParticleMother = pdgSignalParticleArray[iSignPart];
        auto motherIndex = RecoDecay::getMother(mcParticles, mcparticle, pdgParticleMother, true); // check whether mcparticle derives from a particle with pdg = pdgparticlemother, accepting also antiparticle (<- the true parameter)
        if (motherIndex != -1) {
          auto particleMother = mcParticles.rawIteratorAt(motherIndex);
          // just for internal check
          // double mass=particleMother.e()*particleMother.e()-particleMother.pt()*particleMother.pt()-particleMother.pz()*particleMother.pz();
          // filteredTracksMC(mcparticle.pdgCode(),mcparticle.isPhysicalPrimary(),particleMother.pdgCode(),0,motherIndex,0,particleMother.pt(),particleMother.y(),std::sqrt(mass),0);
          if (pdgParticleMother == 310) {
            auto daughtersSlice = mcparticle.template daughters_as<aod::McParticles>();
            int ndaught = daughtersSlice.size(); // might not be accurate in case K0s interact with material before decaying
            if (ndaught != 2)
              ndaught *= -1;
            filteredTracksMC(mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), particleMother.pdgCode(), 0, motherIndex, ndaught, particleMother.pt(), particleMother.y(), 0, 0);
            //  std::cout<<"FOUND K0s, MATCHED!  size array "<<ndaught<<std::endl;
            break;
          }

          int ndaught = 0;
          std::vector<int> indxDaughers;
          if (pdgParticleMother == 421) {
            if (RecoDecay::isMatchedMCGen<true, false>(mcParticles, particleMother, pdgParticleMother, pdgDecayDzero, true, nullptr, 3, &indxDaughers)) {
              ndaught = 2;
              // std::cout<<"########       FOUND D0, MATCHED! pdg: " <<particleMother.pdgCode()<<"################ size array "<<indxDaughers.size()<<std::endl;
            } else
              ndaught = -indxDaughers.size();
          } else if (pdgParticleMother == 4122) {
            if (RecoDecay::isMatchedMCGen<true, false>(mcParticles, particleMother, pdgParticleMother, pdgDecayLc, true, nullptr, 3, &indxDaughers)) {
              ndaught = 3;
            } else
              ndaught = -indxDaughers.size();
          }
          // now check whether the charm hadron is prompt or comes from beauty decay
          std::vector<int> idxBhadMothers;
          if (RecoDecay::getCharmHadronOrigin(mcParticles, particleMother, false, &idxBhadMothers) == RecoDecay::OriginType::NonPrompt) {
            if (idxBhadMothers.size() > 1) {
              LOG(info) << "more than 1 B mother hadron found, should not be: ";
              for (unsigned long iBhM = 0; iBhM < idxBhadMothers.size(); iBhM++) {
                auto particleBhadr = mcParticles.rawIteratorAt(idxBhadMothers[iBhM]);
                LOG(info) << particleBhadr.pdgCode();
              }
            }
            auto particleBhadr = mcParticles.rawIteratorAt(idxBhadMothers[0]);
            // int pdgBhad=particleBhadr.pdgCode();
            filteredTracksMC(mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), particleMother.pdgCode(), particleBhadr.pdgCode(), motherIndex, ndaught, particleMother.pt(), particleMother.y(), particleBhadr.pt(), particleBhadr.y());
          } else {
            filteredTracksMC(mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), particleMother.pdgCode(), 0, motherIndex, ndaught, particleMother.pt(), particleMother.y(), 0, 0);
          }
          break;
        }
        pdgParticleMother = 0;
      }
      if (pdgParticleMother == 0)
        filteredTracksMC(mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), 0, 0, -1, 0, 0, 0, 0, 0);
      // std::cout<<mcparticle.pdgCode()<<std::endl;
    } else {
      filteredTracksMC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }
  }
  void processData(soa::Filtered<TracksWithSelAndDca> const& tracks)
  {
    if (trackPtSampling == 0) {
      for (auto& track : tracks) {
        fillTableData(track);
      }
    } else {
      for (auto& track : lowPtTracks) {
        fillTableData(track);
      }
      for (auto& track : midPtTracks) {
        fillTableData(track);
      }
      for (auto& track : highPtTracks) {
        fillTableData(track);
      }
    }
  }
  PROCESS_SWITCH(FilterTracks, processData, "process data", true);

  void processMC(soa::Filtered<TracksWithSelAndDcaMc> const& tracks, aod::McParticles const& mcParticles)
  {
    if (trackPtSampling == 0) {
      for (auto& track : tracks) {
        fillTableDataMC(track, mcParticles);
      }
    } else {
      for (auto& track : lowPtTracksMC) {
        fillTableDataMC(track, mcParticles);
      }
      for (auto& track : midPtTracksMC) {
        fillTableDataMC(track, mcParticles);
      }
      for (auto& track : highPtTracksMC) {
        fillTableDataMC(track, mcParticles);
      }
    }

    for (auto& mcpart : mcParticles) { // NOTE THAT OF COURSE IN CASE OF SAMPLING THE GEN TABLE WON'T MATCH THE RECO EVEN CONSIDERING EFFICIENCY
      int pdgCode = mcpart.pdgCode();
      // for(int iSignPart=0;iSignPart<3;iSignPart++){

      std::vector<int> indxDaughers;
      float etamax = 0;
      bool isMatchedToSignal = false;
      if (std::abs(pdgCode) == 310) {
        isMatchedToSignal = RecoDecay::isMatchedMCGen<true, false>(mcParticles, mcpart, 310, pdgDecayKzero, true, nullptr, 1, &indxDaughers);
      }
      if (std::abs(pdgCode) == 421) {
        isMatchedToSignal = RecoDecay::isMatchedMCGen<true, false>(mcParticles, mcpart, 421, pdgDecayDzero, true, nullptr, 3, &indxDaughers);
      } else if (std::abs(pdgCode) == 4122) {
        isMatchedToSignal = RecoDecay::isMatchedMCGen<true, false>(mcParticles, mcpart, 4122, pdgDecayLc, true, nullptr, 3, &indxDaughers);
        // std::cout<<"Lc found, matched to MC? "<<isMatchedToSignal<<std::endl;
        // if(!isMatchedToSignal){
        // auto daughtersLxSlice = mcpart.daughters_as<aod::McParticles>();
        //   int ndaught = daughtersLxSlice.size();
        //   for(auto lcDaught : daughtersLxSlice){
        //     std::cout<<"Lc daught, total daught "<<ndaught<<" pdg: "<<lcDaught.pdgCode()<<std::endl;
        //   }
        // }
      }
      if (isMatchedToSignal) {
        for (auto mcpartdaughtIdx : indxDaughers) {
          auto mcPartDaught = mcParticles.rawIteratorAt(mcpartdaughtIdx);
          double eta = std::abs(mcPartDaught.eta());
          if ((eta) > etamax) {
            etamax = eta;
          }
        }
        if (pdgCode == 310) {
          selectedGenParticles(mcpart.pdgCode(), mcpart.mcCollisionId(), 0, mcpart.pt(), mcpart.y(), etamax, 0, 0);
          continue;
        }
        std::vector<int> idxBhadMothers;
        if (RecoDecay::getCharmHadronOrigin(mcParticles, mcpart, false, &idxBhadMothers) == RecoDecay::OriginType::NonPrompt) {
          if (idxBhadMothers.size() > 1) {
            LOG(info) << "loop on gen particles: more than 1 B mother hadron found, should not be: ";
            for (unsigned long iBhM = 0; iBhM < idxBhadMothers.size(); iBhM++) {
              auto particleBhadr = mcParticles.rawIteratorAt(idxBhadMothers[iBhM]);
              LOG(info) << particleBhadr.pdgCode();
            }
          }
          auto particleBhadr = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          // int pdgBhad=particleBhadr.pdgCode();
          selectedGenParticles(mcpart.pdgCode(), mcpart.mcCollisionId(), particleBhadr.pdgCode(), mcpart.pt(), mcpart.y(), etamax, particleBhadr.pt(), particleBhadr.y());
        } else
          selectedGenParticles(mcpart.pdgCode(), mcpart.mcCollisionId(), 0, mcpart.pt(), mcpart.y(), etamax, 0, 0);
      }
      //
    }
  }
  PROCESS_SWITCH(FilterTracks, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FilterTracks>(cfgc)};
}
