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
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;

struct skimmerDalitzEE {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSnn = -1,
  };

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::DalitzEEs> dalitzees;
  Produces<aod::EMPrimaryTracks> emprimarytracks;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.0, "max. TPC n sigma for pion exclusion"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
    },
  };

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(InitContext const&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, typename TTrack>
  bool checkTrack(TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    if (abs(track.tpcNSigmaEl()) > maxTPCNsigmaEl) {
      return false;
    }
    if (abs(track.tpcNSigmaPi()) < maxTPCNsigmaPi) {
      return false;
    }
    return true;
  }

  template <bool isMC, EM_EEPairType pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  void fillPairTable(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if constexpr (pairtype == EM_EEPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(t1) || !checkTrack<isMC>(t2)) {
          continue;
        }
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcaeexy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcaeez = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

        dalitzees(collision.globalIndex(), fNewLabels[t1.globalIndex()], fNewLabels[t2.globalIndex()], v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, dcaeexy, dcaeez, static_cast<int>(pairtype));
        // LOGF(info, "ULS: collision.globalIndex() = %d, fNewLabels[t1.globalIndex()] = %d, fNewLabels[t2.globalIndex()] = %d", collision.globalIndex(), fNewLabels[t1.globalIndex()], fNewLabels[t2.globalIndex()]);
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(t1) || !checkTrack<isMC>(t2)) {
          continue;
        }
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
        float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
        float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
        float dcaeexy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
        float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
        float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
        float dcaeez = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);
        dalitzees(collision.globalIndex(), fNewLabels[t1.globalIndex()], fNewLabels[t2.globalIndex()], v12.Pt(), v12.Eta(), v12.Phi() > 0 ? v12.Phi() : v12.Phi() + TMath::TwoPi(), v12.M(), phiv, dcaeexy, dcaeez, static_cast<int>(pairtype));
        // LOGF(info, "LS: collision.globalIndex() = %d, fNewLabels[t1.globalIndex()] = %d, fNewLabels[t2.globalIndex()] = %d", collision.globalIndex(), fNewLabels[t1.globalIndex()], fNewLabels[t2.globalIndex()]);
      } // end of pairing loop
    }
  }

  template <bool isMC, typename TTracks>
  void fillTrackTable(TTracks const& tracks)
  {
    for (auto& track : tracks) {
      if (!checkTrack<isMC>(track)) {
        continue;
      }

      if (!(fNewLabels.find(track.globalIndex()) != fNewLabels.end())) {
        fNewLabels[track.globalIndex()] = fCounter;

        emprimarytracks(track.collisionId(), track.globalIndex(), track.sign(),
                        track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(),
                        track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                        track.tpcChi2NCl(), track.tpcInnerParam(),
                        track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                        track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                        track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(), track.signed1Pt(), track.cYY(), track.cZZ());

        fCounter++;
      }
    } // end of track loop
  }

  // ============================ FUNCTION DEFINITIONS ====================================================

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& nabs(o2::aod::track::dcaXY) < dca_xy_max&& nabs(o2::aod::track::dcaZ) < dca_z_max&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its;
  using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  std::map<uint64_t, int> fNewLabels;
  int fCounter = 0;
  void processRec(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      // store track info which belongs to pairs. (i.e. only 1 track per event does not enter pair analysis.)
      int npos = 0, nneg = 0;
      for (auto& ptrack : posTracks_per_coll) {
        if (!checkTrack<false>(ptrack)) {
          continue;
        }
        npos++;
      }
      for (auto& ntrack : negTracks_per_coll) {
        if (!checkTrack<false>(ntrack)) {
          continue;
        }
        nneg++;
      }

      if (npos + nneg < 2) {
        continue;
      }

      fillTrackTable<false>(posTracks_per_coll);
      fillTrackTable<false>(negTracks_per_coll);

      fillPairTable<false, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll);  // ULS
      fillPairTable<false, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
      fillPairTable<false, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--

    } // end of collision loop
    fNewLabels.clear();
    fCounter = 0;
  }
  PROCESS_SWITCH(skimmerDalitzEE, processRec, "process reconstructed info only", true);

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      // auto mcCollision = collision.mcCollision();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      // store track info which belongs to pairs. (i.e. only 1 track per event does not enter pair analysis.)
      int npos = 0, nneg = 0;
      for (auto& ptrack : posTracks_per_coll) {
        if (!checkTrack<true>(ptrack)) {
          continue;
        }
        npos++;
      }
      for (auto& ntrack : negTracks_per_coll) {
        if (!checkTrack<true>(ntrack)) {
          continue;
        }
        nneg++;
      }

      if (npos + nneg < 2) {
        continue;
      }

      fillTrackTable<true>(posTracks_per_coll);
      fillTrackTable<true>(negTracks_per_coll);

      fillPairTable<true, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll);  // ULS
      fillPairTable<true, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
      fillPairTable<true, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
    }                                                                                               // end of collision loop
    fNewLabels.clear();
    fCounter = 0;
  }
  PROCESS_SWITCH(skimmerDalitzEE, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerDalitzEE>(cfgc, TaskName{"skimmer-dalitz-ee"})};
}
