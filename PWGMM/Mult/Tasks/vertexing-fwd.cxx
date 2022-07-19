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
#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using FullBCs = soa::Join<aod::BCs, aod::MatchedBCCollisionsSparse>;
using FullCollision = soa::Join<aod::Collisions, aod::McCollisionLabels>;

namespace o2::aod
{
DECLARE_SOA_TABLE(AmbiguousTracksMFT, "AOD", "AMBIGUOUSTRMFT", //! Table for MFT tracks which are not uniquely associated with a collision
                  o2::soa::Index<>, o2::aod::ambiguous::MFTTrackId, o2::aod::ambiguous::BCIdSlice, o2::soa::Marker<2>);
}

struct vertexingfwd {

  /// Could be TEMPORARY: store the vertex, collision, dca, and ambiguous tracks information
  /// into different std::vector to easily handle them later outside the loops
  std::vector<int> vecCollForAmb;        // vector for collisions associated to an ambiguous track
  std::vector<double> vecDCACollForAmb;  // vector for dca collision associated to an ambiguous track
  std::vector<double> vecAmbTrack;       // vector for ambiguous track quantities
  std::vector<double> vecZposCollForAmb; // vector for z vertex of collisions associated to an ambiguous track

  Configurable<float> maxDCAXY{"maxDCAXY", 3.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA

  HistogramRegistry registry{
    "registry",
    {{"TracksDCAXY", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"TracksDCAX", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
     {"TracksDCAY", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
     {"AmbiguousTracksStatus", "; Status; counts", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
     {"NbCollComp", "; NbCollComp", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
     {"NumberOfContributors", "; N_{tr} for vertexing; counts", {HistType::kTH1F, {{100, 0, 100}}}},
     {"CollisionsMatchIndicesMC", "; Rec. minDCA ambitrack coll.ID; Gen. ambitrack coll.ID", {HistType::kTH2F, {{401, -0.5, 1000.5}, {401, -0.5, 1000.5}}}},
     {"TracksDCAXYBest", "; DCA_{xy}^{best} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     //{"TracksDCAXYOther", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"TracksDCAXYBestFalse", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
     {"EfficiencyZvtx", "; z vertex; Efficiency", {HistType::kTProfile, {{100, -30, 30}}}},
     {"DeltaZvtx", "; #delta z (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
     {"DeltaZvtxBest", "; #delta z = z_{best} - z_{true} (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
     //{"DeltaDCAminNcoll", "; N_{coll}; #delta DCA = DCA_{xy}^{min} - DCA_{xy}^{other} (cm); counts", {HistType::kTH2F, {{10, -0.5, 9.5}, {100, -0.5, 9.5}}}},
     {"CorrectMatch", "; Matching value; counts", {HistType::kTH1F, {{4, -0.5, 3.5}}}}}};

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("CorrectMatch"));
    auto* x1 = hstat->GetXaxis();
    x1->SetBinLabel(1, "Incorrect match");
    x1->SetBinLabel(2, "Correct match");
    x1->SetBinLabel(3, "Total w N_{coll} > 0");
    x1->SetBinLabel(4, "Total w N_{coll} #geq 0");

    auto hstatus = registry.get<TH1>(HIST("AmbiguousTracksStatus"));
    auto* x2 = hstatus->GetXaxis();
    x2->SetBinLabel(1, "MFT tracks ");
    x2->SetBinLabel(2, "MFT ambiguous tracks ");
    x2->SetBinLabel(3, "is from primary (all) ");
    x2->SetBinLabel(4, "is from primary ");
    x2->SetBinLabel(5, "not from primary (all) ");
    x2->SetBinLabel(6, "not from primary ");
  }

  int getIndexBestCollision(std::vector<double> vecOfDCA, int method = 0)
  {
    int indice = 0;
    if (vecOfDCA.size() == 0) {
      return -1;
    }
    if (method == 0) {
      indice = std::distance(vecOfDCA.begin(), std::min_element(vecOfDCA.begin(), vecOfDCA.end()));
    }
    return indice;
  }

  void processNew(aod::AmbiguousMFTTracks const& ambitracks, aod::BCs const& bcs, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& tracks, FullCollision const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {

    int ntracks = tracks.size();
    int nambitracks = ambitracks.size();

    registry.fill(HIST("AmbiguousTracksStatus"), 0.0, ntracks);
    registry.fill(HIST("AmbiguousTracksStatus"), 1.0, nambitracks);

    for (auto& ambitrack : ambitracks) {
      vecCollForAmb.clear();
      vecDCACollForAmb.clear();
      vecAmbTrack.clear();
      vecZposCollForAmb.clear();

      double value = 0.0;    // matching value for collision association to an ambiguous track
      double zVtxMCAmbi = 0; // z vertex associated to the mc collision
      int mcCollAmbiID = -1; // mc value for the collision containing the ambiguous track

      auto track = ambitrack.mfttrack_as<MFTTracksLabeled>(); // Obtain the MFT ambiguous track with the MC labels
      // auto extAmbiTrackid = ambitrack.mfttrackId(); // Global index of the MFT ambiguous track

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      mcCollAmbiID = particle.mcCollisionId();
      zVtxMCAmbi = particle.mcCollision().posZ();

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("AmbiguousTracksStatus"), 2.0);
      } else {
        registry.fill(HIST("AmbiguousTracksStatus"), 4.0);
      }

      // Fill the std::vector for ambiguous tracks with the quantities needed
      vecAmbTrack.push_back(track.x());
      vecAmbTrack.push_back(track.y());
      vecAmbTrack.push_back(track.phi());
      vecAmbTrack.push_back(track.tgl());
      vecAmbTrack.push_back(track.signed1Pt());
      vecAmbTrack.push_back(track.z());
      vecAmbTrack.push_back(track.chi2());

      auto bcambis = ambitrack.bc();

      int collCounter = 0;
      for (auto& collision : collisions) {
        uint64_t mostProbableBC = collision.bc().globalBC();
        //uint64_t meanBC = mostProbableBC - std::lround(collision.collisionTime() / (o2::constants::lhc::LHCBunchSpacingNS / 1000));
        //int deltaBC = std::ceil(collision.collisionTimeRes() / (o2::constants::lhc::LHCBunchSpacingNS / 1000) * 4);

        for (auto& bcambi : bcambis) {

          if (bcambi.globalBC() != mostProbableBC) {
            continue;
          }
          collCounter++;

          SMatrix5 tpars(vecAmbTrack[0], vecAmbTrack[1], vecAmbTrack[2], vecAmbTrack[3], vecAmbTrack[4]);
          //            std::vector<double> v1{extAmbiTrack.cXX(), extAmbiTrack.cXY(), extAmbiTrack.cYY(), extAmbiTrack.cPhiX(), extAmbiTrack.cPhiY(),
          //                                   extAmbiTrack.cPhiPhi(), extAmbiTrack.cTglX(), extAmbiTrack.cTglY(), extAmbiTrack.cTglPhi(), extAmbiTrack.cTglTgl(),
          //                                   extAmbiTrack.c1PtX(), extAmbiTrack.c1PtY(), extAmbiTrack.c1PtPhi(), extAmbiTrack.c1PtTgl(), extAmbiTrack.c1Pt21Pt2()};

          std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
          SMatrix55 tcovs(v1.begin(), v1.end());
          o2::track::TrackParCovFwd pars1{vecAmbTrack[5], tpars, tcovs, vecAmbTrack[6]};

          // o2::track::TrackParCovFwd pars1{extAmbiTrack.z(), tpars, tcovs, chi2};
          pars1.propagateToZlinear(collision.posZ()); // track parameters propagation to the position of the z vertex

          const auto dcaX(pars1.getX() - collision.posX());
          const auto dcaY(pars1.getY() - collision.posY());
          auto dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          registry.fill(HIST("TracksDCAXY"), dcaXY);
          registry.fill(HIST("TracksDCAX"), dcaX);
          registry.fill(HIST("TracksDCAY"), dcaY);
          registry.fill(HIST("NumberOfContributors"), collision.numContrib());

          vecDCACollForAmb.push_back(dcaXY);

          if (!collision.has_mcCollision()) {
            continue;
          }

          int mcCollindex = collision.mcCollision().globalIndex();
          vecCollForAmb.push_back(mcCollindex);

          vecZposCollForAmb.push_back(collision.mcCollision().posZ());

          registry.fill(HIST("DeltaZvtx"), collision.mcCollision().posZ() - zVtxMCAmbi);
          break;
        }
      }

      registry.fill(HIST("NbCollComp"), collCounter);
      registry.fill(HIST("CorrectMatch"), 3.0); // counting for ambiguous track with N collisions >=0

      int indexMinDCA = getIndexBestCollision(vecDCACollForAmb, 0); // obtain min value in the stored vector of DCAs
      int indexMCcoll = -1;
      if (indexMinDCA == -1) {
        continue;
      }
      indexMCcoll = vecCollForAmb[indexMinDCA];

      registry.fill(HIST("CollisionsMatchIndicesMC"), mcCollAmbiID, indexMCcoll);

      if (collCounter == 1) {
        printf("strange ambiguous track of mfttrackId %d\n", ambitrack.mfttrackId());
        if (mcCollAmbiID == indexMCcoll) {
          printf("and this is a correct match for the ambiguous track of mfttrackid %d\n", ambitrack.mfttrackId());
        }
      }

      if (vecDCACollForAmb.size() == 0) { // do not use the vector with no collisions
        if (!particle.isPhysicalPrimary()) {
          registry.fill(HIST("AmbiguousTracksStatus"), 5.0);
        }
        if (particle.isPhysicalPrimary()) {
          registry.fill(HIST("AmbiguousTracksStatus"), 3.0);
        }
        continue;
      }
      registry.fill(HIST("DeltaZvtxBest"), vecZposCollForAmb[indexMinDCA] - zVtxMCAmbi);
      //      for (auto& dca : vecDCACollForAmb) {
      //        if (dca != vecDCACollForAmb[indexMinDCA]) {
      //          registry.fill(HIST("DeltaDCAminNcoll"), vecCollForAmb.size(), std::abs(vecDCACollForAmb[indexMinDCA] - dca));
      //          registry.fill(HIST("TracksDCAXYOther"), dca);
      //        }
      //      }
      if (mcCollAmbiID == indexMCcoll) {
        value = 1.0;
        // LOGF(info, " --> Ambitrack correctly associated to collision, dca= %f", vecDCACollForAmb[indexMinDCA]);
      }
      registry.fill(HIST("TracksDCAXYBest"), vecDCACollForAmb[indexMinDCA]);
      registry.fill(HIST("CorrectMatch"), value);
      registry.fill(HIST("EfficiencyZvtx"), zVtxMCAmbi, value);
      registry.fill(HIST("CorrectMatch"), 2.0); // Counting for amibuous track with N collisions > 0

      if (value == 0.0) {
        registry.fill(HIST("TracksDCAXYBestFalse"), vecDCACollForAmb[indexMinDCA]); // Incorrect association with min DCA
      }

    } // ambitracks loop
  }

  PROCESS_SWITCH(vertexingfwd, processNew, "Process ambiguous track DCA", true);

  void processOld(aod::AmbiguousTracksMFT const& ambitracks, aod::BCs const& bcs, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& tracks, FullCollision const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {

    int ntracks = tracks.size();
    int nambitracks = ambitracks.size();

    registry.fill(HIST("AmbiguousTracksStatus"), 0.0, ntracks);
    registry.fill(HIST("AmbiguousTracksStatus"), 1.0, nambitracks);

    for (auto& ambitrack : ambitracks) {
      vecCollForAmb.clear();
      vecDCACollForAmb.clear();
      vecAmbTrack.clear();
      vecZposCollForAmb.clear();

      double value = 0.0;    // matching value for collision association to an ambiguous track
      double zVtxMCAmbi = 0; // z vertex associated to the mc collision
      int mcCollAmbiID = -1; // mc value for the collision containing the ambiguous track

      auto track = ambitrack.mfttrack_as<MFTTracksLabeled>(); // Obtain the MFT ambiguous track with the MC labels
      // auto extAmbiTrackid = ambitrack.mfttrackId(); // Global index of the MFT ambiguous track

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      mcCollAmbiID = particle.mcCollisionId();
      zVtxMCAmbi = particle.mcCollision().posZ();

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("AmbiguousTracksStatus"), 2.0);
      } else {
        registry.fill(HIST("AmbiguousTracksStatus"), 4.0);
      }

      // Fill the std::vector for ambiguous tracks with the quantities needed
      vecAmbTrack.push_back(track.x());
      vecAmbTrack.push_back(track.y());
      vecAmbTrack.push_back(track.phi());
      vecAmbTrack.push_back(track.tgl());
      vecAmbTrack.push_back(track.signed1Pt());
      vecAmbTrack.push_back(track.z());
      vecAmbTrack.push_back(track.chi2());

      auto bcambis = ambitrack.bc();

      int collCounter = 0;
      for (auto& collision : collisions) {
        uint64_t mostProbableBC = collision.bc().globalBC();
        //uint64_t meanBC = mostProbableBC - std::lround(collision.collisionTime() / (o2::constants::lhc::LHCBunchSpacingNS / 1000));
        //int deltaBC = std::ceil(collision.collisionTimeRes() / (o2::constants::lhc::LHCBunchSpacingNS / 1000) * 4);

        for (auto& bcambi : bcambis) {
          if (bcambi.globalBC() != mostProbableBC) {
            continue;
          }
          collCounter++;

          SMatrix5 tpars(vecAmbTrack[0], vecAmbTrack[1], vecAmbTrack[2], vecAmbTrack[3], vecAmbTrack[4]);
          //            std::vector<double> v1{extAmbiTrack.cXX(), extAmbiTrack.cXY(), extAmbiTrack.cYY(), extAmbiTrack.cPhiX(), extAmbiTrack.cPhiY(),
          //                                   extAmbiTrack.cPhiPhi(), extAmbiTrack.cTglX(), extAmbiTrack.cTglY(), extAmbiTrack.cTglPhi(), extAmbiTrack.cTglTgl(),
          //                                   extAmbiTrack.c1PtX(), extAmbiTrack.c1PtY(), extAmbiTrack.c1PtPhi(), extAmbiTrack.c1PtTgl(), extAmbiTrack.c1Pt21Pt2()};

          std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
          SMatrix55 tcovs(v1.begin(), v1.end());
          o2::track::TrackParCovFwd pars1{vecAmbTrack[5], tpars, tcovs, vecAmbTrack[6]};

          // o2::track::TrackParCovFwd pars1{extAmbiTrack.z(), tpars, tcovs, chi2};
          pars1.propagateToZlinear(collision.posZ()); // track parameters propagation to the position of the z vertex

          const auto dcaX(pars1.getX() - collision.posX());
          const auto dcaY(pars1.getY() - collision.posY());
          auto dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          registry.fill(HIST("TracksDCAXY"), dcaXY);
          registry.fill(HIST("TracksDCAX"), dcaX);
          registry.fill(HIST("TracksDCAY"), dcaY);
          registry.fill(HIST("NumberOfContributors"), collision.numContrib());

          vecDCACollForAmb.push_back(dcaXY);

          if (!collision.has_mcCollision()) {
            continue;
          }

          int mcCollindex = collision.mcCollision().globalIndex();
          vecCollForAmb.push_back(mcCollindex);

          vecZposCollForAmb.push_back(collision.mcCollision().posZ());

          registry.fill(HIST("DeltaZvtx"), collision.mcCollision().posZ() - zVtxMCAmbi);
          break;
        }
      }

      registry.fill(HIST("NbCollComp"), collCounter);

      int indexMinDCA = getIndexBestCollision(vecDCACollForAmb, 0); // obtain min value in the stored vector of DCAs
      int indexMCcoll = -1;
      if (indexMinDCA == -1) {
        continue;
      }
      indexMCcoll = vecCollForAmb[indexMinDCA];

      registry.fill(HIST("CollisionsMatchIndicesMC"), mcCollAmbiID, indexMCcoll);

      registry.fill(HIST("CorrectMatch"), 3.0); // counting for amibuous track with N collisions >=0
      if (vecDCACollForAmb.size() == 0) {       // do not use the vector with no collisions
        if (!particle.isPhysicalPrimary()) {
          registry.fill(HIST("AmbiguousTracksStatus"), 5.0);
        }
        if (particle.isPhysicalPrimary()) {
          registry.fill(HIST("AmbiguousTracksStatus"), 3.0);
        }
        continue;
      }
      registry.fill(HIST("DeltaZvtxBest"), vecZposCollForAmb[indexMinDCA] - zVtxMCAmbi);
      //      for (auto& dca : vecDCACollForAmb) {
      //        if (dca != vecDCACollForAmb[indexMinDCA]) {
      //          registry.fill(HIST("DeltaDCAminNcoll"), vecCollForAmb.size(), std::abs(vecDCACollForAmb[indexMinDCA] - dca));
      //          registry.fill(HIST("TracksDCAXYOther"), dca);
      //        }
      //      }
      if (mcCollAmbiID == indexMCcoll) {
        value = 1.0;
        // LOGF(info, " --> Ambitrack correctly associated to collision, dca= %f", vecDCACollForAmb[indexMinDCA]);
      }
      registry.fill(HIST("TracksDCAXYBest"), vecDCACollForAmb[indexMinDCA]);
      registry.fill(HIST("CorrectMatch"), value);
      registry.fill(HIST("EfficiencyZvtx"), zVtxMCAmbi, value);
      registry.fill(HIST("CorrectMatch"), 2.0); // Counting for amibuous track with N collisions > 0
      if (value == 0.0) {
        registry.fill(HIST("TracksDCAXYBestFalse"), vecDCACollForAmb[indexMinDCA]); // Incorrect association with min DCA
      }

    } // ambitracks loop
  }

  PROCESS_SWITCH(vertexingfwd, processOld, "Process ambiguous track DCA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
