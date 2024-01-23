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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent, aod::EMReducedEventsBz>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMReducedEventIds>;
using MyTrack = MyTracks::iterator;

struct skimmerDalitzEE {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSmm = -1,
  };

  SliceCache cache;
  Preslice<MyTracks> perCol = o2::aod::emprimaryelectron::emreducedeventId;

  SliceCache cache_cefp;
  Preslice<aod::EMPrimaryElectrons> perCol_cefp = o2::aod::emprimaryelectron::collisionId;

  Produces<aod::DalitzEEs> dalitzees;
  Produces<o2::aod::DalitzEEEMReducedEventIds> dalitz_ee_eventid;
  Produces<o2::aod::EMReducedEventsNee> event_nee;

  // Configurables
  Configurable<float> maxMee{"maxMee", 0.5, "max. mee to store ee pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
      {"hNele", "hNele;centrality FT0C;Number of electrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
      {"hNpos", "hNpos;centrality FT0C;Number of positrons", {HistType::kTH2F, {{110, 0, 110}, {101, -0.5f, +100.5f}}}},
    },
  };

  void init(InitContext const&) {}

  template <EM_EEPairType pairtype, bool isCEFP, typename TCollision, typename TTracks1, typename TTracks2>
  int fillPairTable(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    int npair = 0;
    if constexpr (pairtype == EM_EEPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMee) { // don't store
          continue;
        }
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());
        float opangle = getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcaeexy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcaeez = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

        if constexpr (isCEFP) {
          dalitzees(collision.globalIndex(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcaeexy, dcaeez, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        } else { // for analysis
          dalitzees(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcaeexy, dcaeez, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        }

        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMee) { // don't store
          continue;
        }
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());
        float opangle = getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcaeexy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcaeez = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

        if constexpr (isCEFP) {
          dalitzees(collision.globalIndex(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcaeexy, dcaeez, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        } else { // for analysis
          dalitzees(collision.collisionId(), t1.globalIndex(), t2.globalIndex(), v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, opangle, dcaeexy, dcaeez, static_cast<int>(pairtype));
          dalitz_ee_eventid(collision.globalIndex());
        }

        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        npair++;
      } // end of pairing loop
    }
    return npair;
  }

  Partition<MyTracks> posTracks = o2::aod::emprimaryelectron::sign > 0;
  Partition<MyTracks> negTracks = o2::aod::emprimaryelectron::sign < 0;
  void processAnalysis(MyCollisions const& collisions, MyTracks const& tracks)
  {
    for (auto& collision : collisions) {
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emreducedeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emreducedeventId, collision.globalIndex(), cache);
      fRegistry.fill(HIST("hNpos"), collision.centFT0C(), posTracks_per_coll.size());
      fRegistry.fill(HIST("hNele"), collision.centFT0C(), negTracks_per_coll.size());
      // LOGF(info, "collision.centFT0C() = %f, posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", collision.centFT0C() , posTracks_per_coll.size(), negTracks_per_coll.size());

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_EEPairType::kULS, false>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_EEPairType::kLSpp, false>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_EEPairType::kLSmm, false>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nee(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzEE, processAnalysis, "Process dalitz ee for analysis", true);

  Partition<aod::EMPrimaryElectrons> posTracks_cefp = o2::aod::emprimaryelectron::sign > 0;
  Partition<aod::EMPrimaryElectrons> negTracks_cefp = o2::aod::emprimaryelectron::sign < 0;
  void processCEFP(soa::Join<aod::Collisions, aod::EMReducedEventsBz> const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    for (auto& collision : collisions) {
      auto posTracks_per_coll = posTracks_cefp->sliceByCached(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache_cefp);
      auto negTracks_per_coll = negTracks_cefp->sliceByCached(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache_cefp);

      int npair_uls = 0, npair_lspp = 0, npair_lsmm = 0;
      npair_uls = fillPairTable<EM_EEPairType::kULS, true>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        npair_lspp = fillPairTable<EM_EEPairType::kLSpp, true>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        npair_lsmm = fillPairTable<EM_EEPairType::kLSmm, true>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
      event_nee(npair_uls, npair_lspp, npair_lsmm);
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerDalitzEE, processCEFP, "Process dalitz ee for CEFP", false); // for central event filter processing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzEE>(cfgc, TaskName{"skimmer-dalitz-ee"})};
}
