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

/// \file taskMultiplicityEstimatorCorrelation.cxx
/// \brief Task for correlating the multiplicity estimator with generated dN/deta
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace check_mc_pv_contr
{
// V0s
DECLARE_SOA_COLUMN(IdCollGen, idCollGen, int);             //! Generated collision index
DECLARE_SOA_COLUMN(ImpParGen, impParGen, float);           //! Generated impact parameter
DECLARE_SOA_COLUMN(XCollGen, xCollGen, float);             //! Generated x coordinate of the collision
DECLARE_SOA_COLUMN(YCollGen, yCollGen, float);             //! Generated y coordinate of the collision
DECLARE_SOA_COLUMN(ZCollGen, zCollGen, float);             //! Generated z coordinate of the collision
DECLARE_SOA_COLUMN(TimeGen, timeGen, float);               //! Generated collision time
DECLARE_SOA_COLUMN(TimeRec, timeRec, float);               //! Reconstructed collision time
DECLARE_SOA_COLUMN(NCharm, nCharm, int);                   //! Number of charm quarks in the collision
DECLARE_SOA_COLUMN(NCharmFromInj, nCharmFromInj, int);     //! Number of charm quarks from injected events
DECLARE_SOA_COLUMN(NPVContributors, nPVContributors, int); //! Number of contributors to the PV
DECLARE_SOA_COLUMN(Centrality, centrality, int);           //! Centrality FT0C
DECLARE_SOA_COLUMN(XCollRec, xCollRec, float);             //! Reconstructed x coordinate of the collision
DECLARE_SOA_COLUMN(YCollRec, yCollRec, float);             //! Reconstructed y coordinate of the collision
DECLARE_SOA_COLUMN(ZCollRec, zCollRec, float);             //! Reconstructed z coordinate of the collision
DECLARE_SOA_COLUMN(BC, Bc, int);                           //! Bunch crossing
} // namespace check_mc_pv_contr

DECLARE_SOA_TABLE(CheckInj, "AOD", "CHECKINJ", //! Table with PID information
                  check_mc_pv_contr::IdCollGen,
                  check_mc_pv_contr::ImpParGen,
                  check_mc_pv_contr::XCollGen,
                  check_mc_pv_contr::YCollGen,
                  check_mc_pv_contr::ZCollGen,
                  check_mc_pv_contr::TimeGen,
                  check_mc_pv_contr::TimeRec,
                  check_mc_pv_contr::NCharm,
                  check_mc_pv_contr::NCharmFromInj,
                  check_mc_pv_contr::NPVContributors,
                  check_mc_pv_contr::Centrality,
                  check_mc_pv_contr::XCollRec,
                  check_mc_pv_contr::YCollRec,
                  check_mc_pv_contr::ZCollRec,
                  check_mc_pv_contr::BC);
} // namespace o2::aod

struct checkMcPvContr {

  Produces<o2::aod::CheckInj> checkInj;

  using TrackWLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using CollisionWLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentNTPVs>;

  PresliceUnsorted<CollisionWLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  Preslice<TrackWLabels> tracksPerCollision = aod::track::collisionId;

  HistogramRegistry registry{"registry", {}};
  std::shared_ptr<TH1> hCharmPerCollImpPar;

  bool isCharm(int pdg)
  {
    if (std::abs(pdg) / 1000 == 4)
      return true;
    if (std::abs(pdg) / 100 == 4)
      return true;
    return false;
  }

  bool isBeauty(int pdg)
  {
    if (std::abs(pdg) / 1000 == 5)
      return true;
    if (std::abs(pdg) / 100 == 5)
      return true;
    return false;
  }

  void init(InitContext&)
  {
    registry.add("hCharmImpPar", ";Impact parameter (fm);Charm counts", {HistType::kTH1F, {{200, 0, 20}}});
    registry.add("hCollImpPar", ";Impact parameter (fm);Counts", {HistType::kTH1F, {{200, 0, 20}}});
    hCharmPerCollImpPar = registry.add<TH1>("hCharmPerCollImpPar", ";Impact parameter (fm);Charm counts per collision", {HistType::kTH1F, {{200, 0, 20}}});

    registry.add("hDeltaX", ";#DeltaX (cm);Counts", {HistType::kTH1F, {{200, -0.01, 0.01}}});
    registry.add("hDeltaY", ";#DeltaY (cm);Counts", {HistType::kTH1F, {{200, -0.01, 0.01}}});
    registry.add("hDeltaZ", ";#DeltaZ (cm);Counts", {HistType::kTH1F, {{200, -0.01, 0.01}}});
  }

  void process(CollisionWLabels const& collisions,
               TrackWLabels const& tracks,
               aod::McParticles const& mcParticles,
               aod::McCollisions const& mcCollisions)
  {
    assert(isCharm(413));
    int splitColls{0};
    for (const auto& mcColl : mcCollisions) {
      const auto collSlice = collisions.sliceBy(colPerMcCollision, mcColl.globalIndex());
      int64_t idxCollMostPVContrib{-1};
      int nPVContribMost{0};
      // First we find the collision with most PV contributors
      for (const auto& collision : collSlice) {
        if (collision.centFT0C() < 20.f) {
          if (collision.numContrib() > nPVContribMost) {
            nPVContribMost = collision.numContrib();
            idxCollMostPVContrib = collision.globalIndex();
          }
        }
      }

      // Then we fill the histogram with the distances of the collisions
      for (const auto& collision : collSlice) {
        const auto collTracks = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
        std::vector<int> charmIds{};
        int fromSignalEv{0};
        for (const auto& track : collTracks) {
          if (track.has_mcParticle()) {
            auto mcPart = track.mcParticle_as<aod::McParticles>();
            for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
              if (isCharm(mother.pdgCode())) { // charm hadron
                registry.fill(HIST("hDeltaX"), collision.posX() - collisions.rawIteratorAt(idxCollMostPVContrib).posX());
                registry.fill(HIST("hDeltaY"), collision.posY() - collisions.rawIteratorAt(idxCollMostPVContrib).posY());
                registry.fill(HIST("hDeltaZ"), collision.posZ() - collisions.rawIteratorAt(idxCollMostPVContrib).posZ());
                if (std::find(charmIds.begin(), charmIds.end(), mother.globalIndex()) == charmIds.end()) {
                  charmIds.push_back(mother.globalIndex());
                  fromSignalEv += static_cast<int>(!mother.fromBackgroundEvent());
                }
                break;
              }
            }
          }
        }
        checkInj(
          mcColl.globalIndex(), mcColl.impactParameter(), mcColl.posX(), mcColl.posY(), mcColl.posZ(), mcColl.t(), collision.collisionTime(),
          charmIds.size(), fromSignalEv, collision.numContrib(), collision.centFT0C(), collision.posX(), collision.posY(), collision.posZ(),
          collision.bcId());
      }
    }
    for (const auto& mcColl : mcCollisions) {
      registry.fill(HIST("hCollImpPar"), mcColl.impactParameter());
    }
    int count_bkg{0}, count_sgn{0};
    for (const auto& mcPart : mcParticles) {
      if (isCharm(mcPart.pdgCode())) { // charm hadron
        if (!mcPart.fromBackgroundEvent()) {
          auto mcCollision = mcPart.mcCollision_as<aod::McCollisions>();
          registry.fill(HIST("hCharmImpPar"), mcCollision.impactParameter());
          count_sgn++;
        } else {
          count_bkg++;
        }
      }
    }
    hCharmPerCollImpPar->Divide(registry.get<TH1>(HIST("hCharmImpPar")).get(), registry.get<TH1>(HIST("hCollImpPar")).get(), 1, 1);
    LOG(info) << "Number of bkgev particles: " << count_bkg;
    LOG(info) << "Number of sngev particles: " << count_sgn;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<checkMcPvContr>(cfgc)};
}
