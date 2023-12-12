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

/// \file tagAndProbeDmesons.cxx
/// \brief Task for tracking efficiency studies with tag-and-probe using 3-prong D-meson decays
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Reconstruction of 2-prong displaced vertices (very good quality and purity)
/// 1) K∓K± for φ from Ds± or D± → φπ± decays
/// 2) π±π± for D± → K∓π±π± decays
/// 3) K∓π± for D0 from D±* → D0π± decays
struct TagTwoProngDisplacedVertices {

  using TracksWithSelAndDca = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection>;

  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && std::abs(aod::track::dcaXY) > 0.002f; // for the tag, we only consider global tracks with dcaXY > 20 microns
  TracksWithSelAndDcaFiltered = soa::Filtered<TracksWithSelAndDca>;

  Partition<TracksWithSelAndDcaFiltered> positive = aod::track::signed1Pt > 0.;
  Partition<TracksWithSelAndDcaFiltered> negative = aod::track::signed1Pt < 0.;

  SliceCache cache;
  Preslice<TracksWithSelAndDcaFiltered> perCollision = aod::track::collisionId;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> vertexer;

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init("http://alice-ccdb.cern.ch");

    vertexer.setPropagateToPCA(true);
    vertexer.setMaxR(200.);
    vertexer.setMaxDZIni(4.);
    vertexer.setMinParamChange(1.e-3);
    vertexer.setMinRelChi2Change(0.9);
    vertexer.setUseAbsDCA(false);
  }

  template<typename TTrack>
  bool isSelectedInvariantMass(const TTrack& trackFirst,
                               const TTrack& trackSecond
                               std::array<float, 2>& masses,
                               float& invMassMin,
                               float& invMassMax)
  {
    std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
    std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond};
    auto invMass2 = RecoDecay::m2(arrMomentum, masses);
    if (invMass2 > invMassMax*invMassMax || invMass2 < invMassMin*invMassMin) {
      return false;
    }
    return true;
  }

  template<typename TTracks>
  void computeCombinatorialSameCharge(TTracks const& tracks, // pool of tracks
                                      std::array<float, 2>& masses,
                                      float& invMassMin,
                                      float& invMassMax
                                      float& bz);
  {
    for (auto trackFirst = tracks.begin(); trackIdPion1 != tracks.end(); ++trackFirst) {
      for (auto trackSecond = trackFirst + 1; trackSecond != tracks.end(); ++trackSecond) {        

        if(!isSelectedInvariantMass(trackFirst, trackSecond, masses, invMassMin, invMassMax)) {
          continue;
        }

        auto trackParCovFirst = getTrackParCov(trackFirst);
        auto trackParCovSecond = getTrackParCov(trackSecond);

        int nVertices{0};
        try {
          nVertices = fitter.process(trackParCovFirst, trackParCovSecond);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices == 0) {
          continue;
        }
        // add topological cuts here (displacement)
      }
    }
  }

  void processPiPi(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                   aod::BCsWithTimestamps const& bc,
                   TracksWithSelAndDcaFiltered tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    auto groupPositive = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    float invMassMin{0.2f}; 
    float invMassMax{1.4f};
    computeCombinatorialSameCharge(groupPositive, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged}, invMassMin, invMassMax);

    auto groupNegative = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

  }

};

// /// Probe third track reconstruction efficiency with different selections
// struct ProbeThirdTrack {

// }

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TagTwoProngDisplacedVertices>(cfgc));
  return workflow;
}