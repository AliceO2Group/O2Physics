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

/// \file taskMcInjection.cxx
/// \brief Task for checking injected events in Pb-Pb MC productions
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino

#include "Common/DataModel/Centrality.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <cstdlib>
#include <memory>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace check_mc_pv_contr
{
// Collisions
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
DECLARE_SOA_COLUMN(Bc, bc, int);                           //! Bunch crossing
// Tracks
DECLARE_SOA_COLUMN(Vx, vx, float);                    // x coordinate of the track production vertex
DECLARE_SOA_COLUMN(Vy, vy, float);                    // y coordinate of the track production vertex
DECLARE_SOA_COLUMN(Vz, vz, float);                    // z coordinate of the track production vertex
DECLARE_SOA_COLUMN(IsFromSignal, isFromSignal, bool); // Whether the track is from the signal event
} // namespace check_mc_pv_contr

DECLARE_SOA_TABLE(CheckInj, "AOD", "CHECKINJ", //! Table with PID information
                  check_mc_pv_contr::IdCollGen,
                  check_mc_pv_contr::ImpParGen,
                  check_mc_pv_contr::XCollGen,
                  check_mc_pv_contr::YCollGen,
                  check_mc_pv_contr::ZCollGen,
                  check_mc_pv_contr::TimeGen,
                  check_mc_pv_contr::XCollRec,
                  check_mc_pv_contr::YCollRec,
                  check_mc_pv_contr::ZCollRec,
                  check_mc_pv_contr::TimeRec,
                  check_mc_pv_contr::NCharm,
                  check_mc_pv_contr::NCharmFromInj,
                  check_mc_pv_contr::NPVContributors,
                  check_mc_pv_contr::Centrality,
                  check_mc_pv_contr::Bc);

DECLARE_SOA_TABLE(TracksInjection, "AOD", "TRKINJ", //! Table with MC labels for tracks
                  check_mc_pv_contr::IdCollGen,
                  check_mc_pv_contr::Vx,
                  check_mc_pv_contr::Vy,
                  check_mc_pv_contr::Vz,
                  check_mc_pv_contr::IsFromSignal);
} // namespace o2::aod

struct HfTaskMcInjection {

  Produces<o2::aod::CheckInj> checkInj;
  Produces<o2::aod::TracksInjection> tracksInj;

  Configurable<double> centMaxForCollDelta{"centMaxForCollDelta", 20., "max. cent. for gen-rec collision position histograms"};
  Configurable<int> nPvContribMaxForCollDelta{"nPvContribMaxForCollDelta", 2000, "max. PV contrib. for gen-rec collision position histograms"};

  std::shared_ptr<TH1> hCharmPerCollImpPar, hCollisions;

  using TrackWLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using CollisionWLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentNTPVs>;

  Preslice<TrackWLabels> tracksPerCollision = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CollisionWLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    AxisSpec impParBins = {200, 0, 20};
    AxisSpec deltaXYbins = {200, -0.05, 0.05};
    AxisSpec deltaZbins = {200, -10, 10};

    registry.add("hCharmImpPar", ";Impact parameter (fm);Charm counts", {HistType::kTH1F, {impParBins}});
    registry.add("hCollImpPar", ";Impact parameter (fm);Counts", {HistType::kTH1F, {impParBins}});
    hCharmPerCollImpPar = registry.add<TH1>("hCharmPerCollImpPar", ";Impact parameter (fm);Charm counts per collision", {HistType::kTH1F, {impParBins}});

    registry.add("hDeltaX", ";#DeltaX (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaY", ";#DeltaY (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaZ", ";#DeltaZ (cm);Counts", {HistType::kTH1F, {{deltaZbins}}});

    registry.add("hDeltaX_NPV_lt", ";#DeltaX (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaY_NPV_lt", ";#DeltaY (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaZ_NPV_lt", ";#DeltaZ (cm);Counts", {HistType::kTH1F, {{deltaZbins}}});

    registry.add("hDeltaX_NPV_gt", ";#DeltaX (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaY_NPV_gt", ";#DeltaY (cm);Counts", {HistType::kTH1F, {{deltaXYbins}}});
    registry.add("hDeltaZ_NPV_gt", ";#DeltaZ (cm);Counts", {HistType::kTH1F, {{deltaZbins}}});

    registry.add("hDeltaXSngBkg", ";#DeltaX (signal/bkg) (cm);Counts", {HistType::kTH1F, {{200, -10, 10}}});
    registry.add("hDeltaYSngBkg", ";#DeltaY (signal/bkg) (cm);Counts", {HistType::kTH1F, {{200, -10, 10}}});
    registry.add("hDeltaZSngBkg", ";#DeltaZ (signal/bkg) (cm);Counts", {HistType::kTH1F, {{200, -20, 20}}});

    hCollisions = registry.add<TH1>("hCollisions", ";;Counts", {HistType::kTH1F, {{2, 0.5, 2.5}}});
    hCollisions->GetXaxis()->SetBinLabel(1, "Generated");
    hCollisions->GetXaxis()->SetBinLabel(2, "Reconstructed");
  }

  bool isCharm(int pdg)
  {
    if (std::abs(pdg) / 1000 == PDG_t::kCharm) // o2-linter: disable=magic-number (get thousands digit)
      return true;
    if (std::abs(pdg) / 100 == PDG_t::kCharm) // o2-linter: disable=magic-number (get hundreds digit)
      return true;
    return false;
  }

  bool isBeauty(int pdg) // if needed in the future
  {
    if (std::abs(pdg) / 1000 == PDG_t::kBottom) // o2-linter: disable=magic-number (get thousands digit)
      return true;
    if (std::abs(pdg) / 100 == PDG_t::kBottom) // o2-linter: disable=magic-number (get hundreds digit)
      return true;
    return false;
  }

  void process(CollisionWLabels const& collisions,
               TrackWLabels const& tracks,
               aod::McParticles const& mcParticles,
               aod::McCollisions const& mcCollisions)
  {
    for (const auto& mcColl : mcCollisions) {
      registry.fill(HIST("hCollImpPar"), mcColl.impactParameter());
      const auto mcPartColl = mcParticles.sliceBy(perMcCollision, mcColl.globalIndex());
      double xAvgSgn{0.}, yAvgSgn{0.}, zAvgSgn{0.};
      double xAvgBkg{0.}, yAvgBkg{0.}, zAvgBkg{0.};
      int nSgn{0}, nBkg{0};
      for (const auto& mcPart : mcPartColl) {
        if (isCharm(mcPart.pdgCode())) { // charm hadron
          registry.fill(HIST("hCharmImpPar"), mcColl.impactParameter());
        }
        if (mcPart.fromBackgroundEvent()) {
          xAvgBkg += mcPart.vx();
          yAvgBkg += mcPart.vy();
          zAvgBkg += mcPart.vz();
          nBkg++;
          tracksInj(mcPart.mcCollisionId(), mcPart.vx(), mcPart.vy(), mcPart.vz(), false);
        } else {
          xAvgSgn += mcPart.vx();
          yAvgSgn += mcPart.vy();
          zAvgSgn += mcPart.vz();
          nSgn++;
          tracksInj(mcPart.mcCollisionId(), mcPart.vx(), mcPart.vy(), mcPart.vz(), true);
        }
      }
      registry.fill(HIST("hDeltaXSngBkg"), xAvgSgn / nSgn - xAvgBkg / nBkg);
      registry.fill(HIST("hDeltaYSngBkg"), yAvgSgn / nSgn - yAvgBkg / nBkg);
      registry.fill(HIST("hDeltaZSngBkg"), zAvgSgn / nSgn - zAvgBkg / nBkg);

      const auto collSlice = collisions.sliceBy(colPerMcCollision, mcColl.globalIndex());

      // Then we fill the histogram with the distances of the collisions
      for (const auto& collision : collSlice) {
        const auto collTracks = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
        int fromSignalEv{0};
        if (collision.centFT0C() < centMaxForCollDelta) {
          registry.fill(HIST("hDeltaX"), collision.posX() - collision.mcCollision().posX());
          registry.fill(HIST("hDeltaY"), collision.posY() - collision.mcCollision().posY());
          registry.fill(HIST("hDeltaZ"), collision.posZ() - collision.mcCollision().posZ());

          if (collision.numContrib() > nPvContribMaxForCollDelta) {
            registry.fill(HIST("hDeltaX_NPV_gt"), collision.posX() - collision.mcCollision().posX());
            registry.fill(HIST("hDeltaY_NPV_gt"), collision.posY() - collision.mcCollision().posY());
            registry.fill(HIST("hDeltaZ_NPV_gt"), collision.posZ() - collision.mcCollision().posZ());
          } else {
            registry.fill(HIST("hDeltaX_NPV_lt"), collision.posX() - collision.mcCollision().posX());
            registry.fill(HIST("hDeltaY_NPV_lt"), collision.posY() - collision.mcCollision().posY());
            registry.fill(HIST("hDeltaZ_NPV_lt"), collision.posZ() - collision.mcCollision().posZ());
          }
        }
        std::unordered_set<int> charmIds{};
        for (const auto& track : collTracks) {
          if (track.has_mcParticle()) {
            auto mcPart = track.mcParticle_as<aod::McParticles>();
            for (const auto& mother : mcPart.mothers_as<aod::McParticles>()) {
              if (isCharm(mother.pdgCode())) { // charm hadron
                if (!charmIds.contains(mother.globalIndex())) {
                  charmIds.emplace(mother.globalIndex());
                  fromSignalEv += static_cast<int>(!mother.fromBackgroundEvent());
                }
                break;
              }
            }
          }
        }
        checkInj(
          mcColl.globalIndex(), mcColl.impactParameter(),
          mcColl.posX(), mcColl.posY(), mcColl.posZ(), mcColl.t(),
          collision.posX(), collision.posY(), collision.posZ(), collision.collisionTime(),
          charmIds.size(), fromSignalEv, collision.numContrib(), collision.centFT0C(), collision.bcId());
      }
    }

    hCharmPerCollImpPar->Divide(registry.get<TH1>(HIST("hCharmImpPar")).get(), registry.get<TH1>(HIST("hCollImpPar")).get(), 1, 1);
    hCollisions->Fill(1, mcCollisions.size());
    hCollisions->Fill(2, collisions.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcInjection>(cfgc)};
}
