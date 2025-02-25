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
///
/// EMCAL cell track matching task
///
/// \file emcalCellTrackMatcher.cxx
/// \brief Task that matches EMCal cells and tracks
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) Goethe-University
///

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>
#include <random>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "DetectorsBase/GeometryManager.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/CellLabel.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "EMCALBase/Geometry.h"
#include "EMCALBase/ClusterFactory.h"
#include "EMCALBase/NonlinearityHandler.h"
#include "EMCALReconstruction/Clusterizer.h"
#include "PWGJE/Core/EmcalMatchingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyGlobTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TracksDCA, o2::aod::TracksCov>;
using BcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredTracks = o2::soa::Filtered<MyGlobTracks>;
using McCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using FilteredCells = o2::soa::Filtered<aod::Calos>;

struct EmcalCellTrackMatcher {
  Produces<o2::aod::EMCALCellTracks> cellMatchedTracks;

  // Options for the clusterization
  // 1 corresponds to EMCAL cells based on the Run2 definition.
  Configurable<int> selectedCellType{"selectedCellType", 1, "EMCAL Cell type"};
  Configurable<float> trackMinPt{"trackMinPt", 0.3f, "Minimum pT for tracks to perform track matching, to reduce computing time. Tracks below a certain pT will be loopers anyway."};
  Configurable<float> trackDCAz{"trackDCAz", 5.f, "Maximum DCAz of a track to to be considered primary and perform the track matching."};
  Configurable<float> trackDCAxy{"trackDCAxy", 5.f, "Maximum DCAxy of a track to to be considered primary and perform the track matching."};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.01011162697, "Max matching distance track-cell. Default value is equal to half a cell diagonal."};
  Configurable<int> maxNumberTracks{"maxNumberTracks", 1000, "Number of tracks we expect per BC that point to EMCal. This is for optimization."};
  Configurable<bool> useAlignmentFromCCDB{"useAlignmentFromCCDB", true, "States if alignment objects should be used from CCDB"};

  // Require EMCAL cells (CALO type 1)
  Partition<aod::Calos> emcalCells = aod::calo::caloType == selectedCellType;
  Filter emccellfilter = aod::calo::caloType == selectedCellType;

  // Filter for the tracks
  const float trackNotOnEMCal = -900;
  Filter trackFilter = (aod::track::pt >= trackMinPt) && (aod::track::trackEtaEmcal > trackNotOnEMCal) && (aod::track::trackPhiEmcal > trackNotOnEMCal) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (nabs(aod::track::dcaXY) < trackDCAxy) && (nabs(aod::track::dcaZ) < trackDCAz);

  // CDB service (for geometry)
  Service<o2::ccdb::BasicCCDBManager> mCcdbManager;

  // QA
  o2::framework::HistogramRegistry mHistManager{"EmcalCellTrackMatcherQAHistograms"};

  // EMCal geometry
  o2::emcal::Geometry* geometry;
  const int16_t nCellsInEMCal = 12288;
  int nCellsTotal;
  std::vector<std::vector<float>> vecLookUpCellPosition;

  // Preslices
  Preslice<MyGlobTracks> perCollision = o2::aod::track::collisionId;
  PresliceUnsorted<CollEventSels> collisionsPerFoundBC = aod::evsel::foundBCId;
  Preslice<aod::Calos> cellsPerFoundBC = aod::calo::bcId;
  Preslice<aod::SortedTracks> tracksPerFoundBC = aod::emcalcluster::bcId;

  void init(InitContext const&)
  {
    LOG(debug) << "Start init!";
    // NOTE: The geometry manager isn't necessary just to load the EMCAL geometry.
    //       However, it _is_ necessary for loading the misalignment matrices as of September 2020
    //       Eventually, those matrices will be moved to the CCDB, but it's not yet ready.

    // The geomatrices from the CCDB are needed in order to calculate the cluster kinematics
    const char* ccdburl = "http://alice-ccdb.cern.ch"; // Possibly the parameter might become configurable, in order to be able to load new EMCAL calibration objects
    mCcdbManager->setURL(ccdburl);
    mCcdbManager->setCaching(true);
    mCcdbManager->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      mCcdbManager->get<TGeoManager>("GLO/Config/Geometry");
    }
    LOG(debug) << "After load geometry!";
    const int runNumberForGeom = 223409;
    geometry = o2::emcal::Geometry::GetInstanceFromRunNumber(runNumberForGeom);
    nCellsTotal = geometry->GetNCells();
    LOG(debug) << "Completed init!";

    // Setup QA hists.
    using O2HistType = o2::framework::HistType;
    o2::framework::AxisSpec etaAxis{1000, -0.8, 0.8, "#it{#eta}"},
      phiAxis{1000, 0, 2 * 3.14159, "#it{#varphi} (rad)"},
      dEtaAxis{100, -0.015, 0.015, "#Delta#it{#eta}"},
      dPhiAxis{100, -0.015, 0.015, "#Delta#it{#varphi} (rad)"},
      cellAxis{nCellsTotal, -0.5, nCellsTotal - 0.5, "cellID"},
      etaIndexAxis{96, -0.5, 95.5, "#it{#eta} ID"},
      phiIndexAxis{208, -0.5, 207.5, "#it{#varphi} ID"},
      smAxis{20, -0.5, 19.5, "SM"};

    mHistManager.add("hvTrackEtaPhiEMCal", "hvTrackEtaPhiEMCal", O2HistType::kTH2D, {etaAxis, phiAxis});
    mHistManager.add("hTrackCellDiffEMCal", "hTrackCellDiffEMCal", O2HistType::kTH2D, {dEtaAxis, dPhiAxis});
    mHistManager.add("hTrackCellDiffDCal", "hTrackCellDiffDCal", O2HistType::kTH2D, {dEtaAxis, dPhiAxis});
    mHistManager.add("hNTracksPerCell", "hNTracksPerCell", O2HistType::kTH2D, {etaIndexAxis, phiIndexAxis});
    mHistManager.add("hEMCalEtaPhi3DMap", "hEMCalEtaPhi3DMap", O2HistType::kTH2D, {etaIndexAxis, phiIndexAxis});
    mHistManager.add("pEMCalRadiusSM", "pEMCalRadiusSM", O2HistType::kTProfile, {smAxis});
    mHistManager.add("pEMCalRadiusEta", "pEMCalRadiusEta", O2HistType::kTProfile, {etaIndexAxis});
    mHistManager.add("pEMCalRadiusPhi", "pEMCalRadiusPhi", O2HistType::kTProfile, {phiIndexAxis});
    const int nSM = 20; // we have 20 SM
    for (int ism = 0; ism < nSM; ++ism) {
      mHistManager.add(Form("SMwise/h2TrackCellDiffSM%d", ism), Form("#Delta#it{#eta} and #Delta#it{#varphi} between tracks and cells for SM %d", ism), O2HistType::kTH2D, {dEtaAxis, dPhiAxis});
    }

    if (useAlignmentFromCCDB.value) {
      geometry->SetMisalMatrixFromCcdb("Users/m/mhemmer/EMCAL/Config/GeometryAligned", 10000);
    }
    double dist = 0; // we will check the position on the surface of the cell
    double lxyzi[3] = {0., 0., 0.}, xyzi[3] = {0., 0., 0.};
    vecLookUpCellPosition = std::vector<std::vector<float>>(nCellsTotal, std::vector<float>(2, -999.f));
    for (int iCell = 0; iCell < nCellsTotal; ++iCell) {
      // get the local coordinates of the cell
      try {
        geometry->RelPosCellInSModule(iCell, dist).GetCoordinates(lxyzi[0], lxyzi[1], lxyzi[2]);
      } catch (o2::emcal::InvalidCellIDException& e) {
        LOG(error) << e.what();
        continue;
      }
      // Now get the global coordinate
      int SM = geometry->GetSuperModuleNumber(iCell);
      geometry->GetGlobal(lxyzi, xyzi, SM);
      auto [iRow, iCol] = geometry->GlobalRowColFromIndex(iCell);
      math_utils::Point3D<float> pos(xyzi[0], xyzi[1], xyzi[2]);
      vecLookUpCellPosition[iCell][0] = pos.Eta();
      vecLookUpCellPosition[iCell][1] = RecoDecay::constrainAngle(pos.Phi());
      mHistManager.fill(HIST("hEMCalEtaPhi3DMap"), iCol, iRow, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusSM"), SM, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusEta"), iCol, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusPhi"), iRow, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
    }
  }

  template <int supermoduleID>
  void supermoduleHistHelper(double dEta, double dPhi)
  {
    static constexpr std::string_view CellTrackDiffHistSM[20] = {"SMwise/h2TrackCellDiffSM0", "SMwise/h2TrackCellDiffSM1", "SMwise/h2TrackCellDiffSM2", "SMwise/h2TrackCellDiffSM3", "SMwise/h2TrackCellDiffSM4", "SMwise/h2TrackCellDiffSM5", "SMwise/h2TrackCellDiffSM6", "SMwise/h2TrackCellDiffSM7", "SMwise/h2TrackCellDiffSM8", "SMwise/h2TrackCellDiffSM9", "SMwise/h2TrackCellDiffSM10", "SMwise/h2TrackCellDiffSM11", "SMwise/h2TrackCellDiffSM12", "SMwise/h2TrackCellDiffSM13", "SMwise/h2TrackCellDiffSM14", "SMwise/h2TrackCellDiffSM15", "SMwise/h2TrackCellDiffSM16", "SMwise/h2TrackCellDiffSM17", "SMwise/h2TrackCellDiffSM18", "SMwise/h2TrackCellDiffSM19"};
    mHistManager.fill(HIST(CellTrackDiffHistSM[supermoduleID]), dEta, dPhi);
  }

  void fillSupermoduleHistograms(int supermoduleID, double dEta, double dPhi)
  {
    // workaround to have the histogram names per supermodule
    // handled as constexpr
    switch (supermoduleID) {
      case 0:
        supermoduleHistHelper<0>(dEta, dPhi);
        break;
      case 1:
        supermoduleHistHelper<1>(dEta, dPhi);
        break;
      case 2:
        supermoduleHistHelper<2>(dEta, dPhi);
        break;
      case 3:
        supermoduleHistHelper<3>(dEta, dPhi);
        break;
      case 4:
        supermoduleHistHelper<4>(dEta, dPhi);
        break;
      case 5:
        supermoduleHistHelper<5>(dEta, dPhi);
        break;
      case 6:
        supermoduleHistHelper<6>(dEta, dPhi);
        break;
      case 7:
        supermoduleHistHelper<7>(dEta, dPhi);
        break;
      case 8:
        supermoduleHistHelper<8>(dEta, dPhi);
        break;
      case 9:
        supermoduleHistHelper<9>(dEta, dPhi);
        break;
      case 10:
        supermoduleHistHelper<10>(dEta, dPhi);
        break;
      case 11:
        supermoduleHistHelper<11>(dEta, dPhi);
        break;
      case 12:
        supermoduleHistHelper<12>(dEta, dPhi);
        break;
      case 13:
        supermoduleHistHelper<13>(dEta, dPhi);
        break;
      case 14:
        supermoduleHistHelper<14>(dEta, dPhi);
        break;
      case 15:
        supermoduleHistHelper<15>(dEta, dPhi);
        break;
      case 16:
        supermoduleHistHelper<16>(dEta, dPhi);
        break;
      case 17:
        supermoduleHistHelper<17>(dEta, dPhi);
        break;
      case 18:
        supermoduleHistHelper<18>(dEta, dPhi);
        break;
      case 19:
        supermoduleHistHelper<19>(dEta, dPhi);
        break;
      default:
        LOG(info) << "Case Supermodule " << supermoduleID;
        break;
    }
  }

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processFull(BcEvSels const& bcs, CollEventSels const& collisions, FilteredTracks const& tracks, aod::Calos const& cells)
  {
    LOG(debug) << "Starting process full.";
    int previousCollisionId = 0; // Collision ID of the last unique BC. Needed to skip unordered collisions to ensure ordered collisionIds in the cluster table
    int nBCsProcessed = 0;
    int nCellsProcessed = 0;

    std::vector<int> matchIndexCell;
    matchIndexCell.reserve(maxNumberTracks); // reserve enough space for better performance

    // outer vector is size as number of ALL calocells in current DF, inside vector will store trackIDs
    std::vector<std::vector<int>> vecCellToTracks(cells.size(), std::vector<int>());
    for (const auto& bc : bcs) {
      LOG(debug) << "Next BC";

      // Get the collisions matched to the BC using foundBCId of the collision
      auto collisionsInFoundBC = collisions.sliceBy(collisionsPerFoundBC, bc.globalIndex());
      auto cellsInBC = emcalCells.sliceBy(cellsPerFoundBC, bc.globalIndex());

      if (!cellsInBC.size()) {
        LOG(debug) << "No cells found for BC";
        continue;
      }
      // We need bcs which only have a single collision matched to them,
      // so that we can get the tracks from that collision
      if (collisionsInFoundBC.size() == 1) {
        // For the matching of a track to a cell we will use a KDTree approach from emcmatchingutilities
        // where cells are matched to tracks. For this we need vectors of cell and track eta and phi.
        std::vector<float> vTrackPhi;
        std::vector<float> vTrackEta;
        std::vector<int> vTrackIndex;
        std::vector<float> vCellPhi;
        std::vector<float> vCellEta;
        std::vector<int> vCellIndex;
        std::vector<int16_t> vTowerIndex;
        LOG(detail) << "Number of cells for BC (CF): " << cellsInBC.size();
        nCellsProcessed += cellsInBC.size();
        // dummy loop to get the first collision
        for (const auto& col : collisionsInFoundBC) {
          if (previousCollisionId > col.globalIndex()) {
            continue;
          }
          previousCollisionId = col.globalIndex();
          if (col.foundBCId() == bc.globalIndex()) {
            // loop over cells and fill cell information
            vCellPhi.reserve(cellsInBC.size());
            vCellEta.reserve(cellsInBC.size());
            vCellIndex.reserve(cellsInBC.size());
            vTowerIndex.reserve(cellsInBC.size());
            for (const auto& cell : cellsInBC) {
              int16_t cellNumber = cell.cellNumber();
              if (cellNumber < 0 || cellNumber >= nCellsTotal) {
                continue;
              }
              vCellEta.emplace_back(vecLookUpCellPosition[cellNumber][0]);
              vCellPhi.emplace_back(vecLookUpCellPosition[cellNumber][1]);
              vCellIndex.emplace_back(cell.globalIndex());
              vTowerIndex.emplace_back(cellNumber);
            }

            // loop over tracks and fill track info
            auto groupedTracks = tracks.sliceBy(perCollision, col.globalIndex());
            vTrackPhi.reserve(groupedTracks.size());
            vTrackEta.reserve(groupedTracks.size());
            vTrackIndex.reserve(groupedTracks.size());
            for (const auto& track : groupedTracks) {
              if (!track.isGlobalTrack()) {
                LOG(info) << "Track is not global and was not filtered!";
                continue;
              }
              vTrackPhi.emplace_back(RecoDecay::constrainAngle(track.trackPhiEmcal()));
              vTrackEta.emplace_back(track.trackEtaEmcal());
              vTrackIndex.emplace_back(track.globalIndex());
            } // end of loop over track

            matchIndexCell.clear();
            emcmatchingutilities::matchCellsAndTracks(vCellPhi, vCellEta, vTrackPhi, vTrackEta, maxMatchingDistance, matchIndexCell);
            for (size_t iTrack = 0; iTrack < matchIndexCell.size(); iTrack++) {
              int cellIDLocal = matchIndexCell[iTrack];
              if (cellIDLocal != -1) {
                int trackID = vTrackIndex[iTrack];
                int cellID = vCellIndex[cellIDLocal];
                vecCellToTracks[cellID].emplace_back(trackID);
                mHistManager.fill(HIST("hvTrackEtaPhiEMCal"), vTrackEta[iTrack], vTrackPhi[iTrack]);
                double dEta = vTrackEta[iTrack] - vCellEta[cellIDLocal];
                double dPhi = vTrackPhi[iTrack] - vCellPhi[cellIDLocal];
                if (nCellsInEMCal > vTowerIndex[cellIDLocal]) {
                  mHistManager.fill(HIST("hTrackCellDiffEMCal"), dEta, dPhi);
                } else {
                  mHistManager.fill(HIST("hTrackCellDiffDCal"), dEta, dPhi);
                }
                auto [iRow, iCol] = geometry->GlobalRowColFromIndex(vTowerIndex[cellIDLocal]);
                mHistManager.fill(HIST("hNTracksPerCell"), iCol, iRow);
                fillSupermoduleHistograms(geometry->GetSuperModuleNumber(vTowerIndex[cellIDLocal]), dEta, dPhi);
              }
            }
          }
        } // end of loop over current collision mathed to the BC
      } // if (collisionsInFoundBC.size() == 1) {
      LOG(debug) << "Done with process BC.";
      nBCsProcessed++;
    } // end of bc loop
    for (const auto& trackIDs : vecCellToTracks) {
      cellMatchedTracks(trackIDs);
    }
    LOG(detail) << "Processed " << nBCsProcessed << " BCs with " << nCellsProcessed << " cells";
  }
  PROCESS_SWITCH(EmcalCellTrackMatcher, processFull, "run full analysis", true);

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processSorted(o2::aod::SortedTracks const& tracks, aod::Calos const& cells)
  {
    LOG(debug) << "Starting process sorted tracks.";

    // outer vector is size as number of ALL calo cells in current DF, inside vector will store trackIDs
    std::vector<std::vector<int>> vecCellToTracks(cells.size(), std::vector<int>());

    auto track = tracks.begin();
    auto cell = emcalCells.begin();

    const auto cellEnd = emcalCells.end();
    const auto trackEnd = tracks.end();

    std::vector<int> matchIndexCell;
    matchIndexCell.reserve(maxNumberTracks); // reserve enough space for better performance

    // For the matching of a track to a cell we will use a KDTree approach from emcmatchingutilities
    // where cells are matched to tracks. For this we need vectors of cell and track eta and phi.
    std::vector<float> vTrackPhi;
    std::vector<float> vTrackEta;
    std::vector<int> vTrackIndex;
    std::vector<float> vCellPhi;
    std::vector<float> vCellEta;
    std::vector<int> vCellIndex;
    std::vector<int16_t> vTowerIndex;
    std::vector<int> vSMIndex;
    vTrackPhi.reserve(maxNumberTracks);
    vTrackEta.reserve(maxNumberTracks);
    vTrackIndex.reserve(maxNumberTracks);
    vCellPhi.reserve(nCellsTotal);
    vCellEta.reserve(nCellsTotal);
    vCellIndex.reserve(nCellsTotal);
    vTowerIndex.reserve(nCellsTotal);
    vSMIndex.reserve(nCellsTotal);
    while (cell != cellEnd && track != trackEnd) {
      auto cellBC = cell.bcId();
      auto trackBC = track.bcId();

      if (cellBC == trackBC) {
        auto currentBCID = cellBC;

        while (cell != cellEnd && cellBC == currentBCID) {
          int16_t cellNumber = cell.cellNumber();
          if (cellNumber < 0 || cellNumber >= nCellsTotal) {
            LOG(info) << "cell number " << cellNumber << " not within EMCal!";
            continue;
          }
          vCellEta.emplace_back(vecLookUpCellPosition[cellNumber][0]);
          vCellPhi.emplace_back(vecLookUpCellPosition[cellNumber][1]);
          vCellIndex.emplace_back(cell.globalIndex());
          vTowerIndex.emplace_back(cellNumber);
          vSMIndex.emplace_back(geometry->GetSuperModuleNumber(cellNumber));
          if (++cell != cellEnd) {
            cellBC = cell.bcId();
          }
        }

        while (track != trackEnd && trackBC == currentBCID) {
          vTrackPhi.emplace_back(track.phi());
          vTrackEta.emplace_back(track.eta());
          vTrackIndex.emplace_back(track.trackId());
          if (++track != trackEnd) {
            trackBC = track.bcId();
          }
        }
        matchIndexCell.clear();
        emcmatchingutilities::matchCellsAndTracks(vCellPhi, vCellEta, vTrackPhi, vTrackEta, maxMatchingDistance, matchIndexCell);
        for (size_t iTrack = 0; iTrack < matchIndexCell.size(); ++iTrack) {
          float trackPhi = vTrackPhi[iTrack];
          float trackEta = vTrackEta[iTrack];
          int cellIDLocal = matchIndexCell[iTrack];
          if (cellIDLocal != -1) {
            int trackID = vTrackIndex[iTrack];
            int cellID = vCellIndex[cellIDLocal];
            vecCellToTracks[cellID].emplace_back(trackID);
            mHistManager.fill(HIST("hvTrackEtaPhiEMCal"), trackEta, trackPhi);
            double dEta = trackEta - vCellEta[cellIDLocal];
            double dPhi = trackPhi - vCellPhi[cellIDLocal];
            if (nCellsInEMCal > vTowerIndex[cellIDLocal]) {
              mHistManager.fill(HIST("hTrackCellDiffEMCal"), dEta, dPhi);
            } else {
              mHistManager.fill(HIST("hTrackCellDiffDCal"), dEta, dPhi);
            }
            auto [iRow, iCol] = geometry->GlobalRowColFromIndex(vTowerIndex[cellIDLocal]);
            mHistManager.fill(HIST("hNTracksPerCell"), iCol, iRow);
            fillSupermoduleHistograms(vSMIndex[cellIDLocal], dEta, dPhi);
          }
        }
        vTrackPhi.clear();
        vTrackEta.clear();
        vTrackIndex.clear();
        vCellPhi.clear();
        vCellEta.clear();
        vCellIndex.clear();
        vTowerIndex.clear();
        vSMIndex.clear();
      } else if (cellBC < trackBC) {
        while (cell != cellEnd && cellBC < trackBC) {
          if (++cell != cellEnd) {
            cellBC = cell.bcId();
          }
        }
      } else {
        while (track != trackEnd && trackBC < cellBC) {
          if (++track != trackEnd) {
            trackBC = track.bcId();
          }
        }
      }
      if (cell == cellEnd || track == trackEnd) {
        break;
      }
    } // while (cell != cellEnd && track != trackEnd)

    for (const auto& trackIDs : vecCellToTracks) {
      cellMatchedTracks(trackIDs);
    }
  }
  PROCESS_SWITCH(EmcalCellTrackMatcher, processSorted, "run analysis with tracks sorted according to BCId", false);

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processSortedWithSlices(BcEvSels const& bcs, o2::aod::SortedTracks const& tracks, aod::Calos const& cells)
  {
    LOG(debug) << "Starting process sorted with slices.";

    std::vector<int> matchIndexCell;
    matchIndexCell.reserve(maxNumberTracks); // reserve enough space for better performance

    // outer vector is size as number of ALL calocells in current DF, inside vector will store trackIDs
    std::vector<std::vector<int>> vecCellToTracks(cells.size(), std::vector<int>());
    for (const auto& bc : bcs) {

      // Get the collisions matched to the BC using foundBCId of the collision
      auto cellsInBC = emcalCells.sliceBy(cellsPerFoundBC, bc.globalIndex());
      auto tracksInBc = tracks.sliceBy(tracksPerFoundBC, bc.globalIndex());

      if (!cellsInBC.size()) {
        LOG(debug) << "No cells found for current BC";
        continue;
      }
      if (!tracksInBc.size()) {
        LOG(debug) << "No tracks found for current BC";
        continue;
      }
      // For the matching of a track to a cell we will use a KDTree approach from emcmatchingutilities
      // where cells are matched to tracks. For this we need vectors of cell and track eta and phi.
      std::vector<float> vTrackPhi;
      std::vector<float> vTrackEta;
      std::vector<int> vTrackIndex;
      std::vector<float> vCellPhi;
      std::vector<float> vCellEta;
      std::vector<int> vCellIndex;
      std::vector<int16_t> vTowerIndex;

      // loop over cells and fill cell information
      vCellPhi.reserve(cellsInBC.size());
      vCellEta.reserve(cellsInBC.size());
      vCellIndex.reserve(cellsInBC.size());
      vTowerIndex.reserve(cellsInBC.size());
      for (const auto& cell : cellsInBC) {
        int16_t cellNumber = cell.cellNumber();
        vCellEta.emplace_back(vecLookUpCellPosition[cellNumber][0]);
        vCellPhi.emplace_back(vecLookUpCellPosition[cellNumber][1]);
        vCellIndex.emplace_back(cell.globalIndex());
        vTowerIndex.emplace_back(cellNumber);
      }

      vTrackPhi.reserve(tracksInBc.size());
      vTrackEta.reserve(tracksInBc.size());
      vTrackIndex.reserve(tracksInBc.size());
      for (const auto& track : tracksInBc) {
        vTrackPhi.emplace_back(track.phi());
        vTrackEta.emplace_back(track.eta());
        vTrackIndex.emplace_back(track.globalIndex());
      } // end of loop over track

      matchIndexCell.clear();
      emcmatchingutilities::matchCellsAndTracks(vCellPhi, vCellEta, vTrackPhi, vTrackEta, maxMatchingDistance, matchIndexCell);
      for (size_t iTrack = 0; iTrack < matchIndexCell.size(); iTrack++) {
        int cellIDLocal = matchIndexCell[iTrack];
        if (cellIDLocal != -1) {
          int trackID = vTrackIndex[iTrack];
          int cellID = vCellIndex[cellIDLocal];
          vecCellToTracks[cellID].emplace_back(trackID);
          mHistManager.fill(HIST("hvTrackEtaPhiEMCal"), vTrackEta[iTrack], vTrackPhi[iTrack]);
          double dEta = vTrackEta[iTrack] - vCellEta[cellIDLocal];
          double dPhi = vTrackPhi[iTrack] - vCellPhi[cellIDLocal];
          if (nCellsInEMCal > vTowerIndex[cellIDLocal]) {
            mHistManager.fill(HIST("hTrackCellDiffEMCal"), dEta, dPhi);
          } else {
            mHistManager.fill(HIST("hTrackCellDiffDCal"), dEta, dPhi);
          }
          auto [iRow, iCol] = geometry->GlobalRowColFromIndex(vTowerIndex[cellIDLocal]);
          mHistManager.fill(HIST("hNTracksPerCell"), iCol, iRow);
          fillSupermoduleHistograms(geometry->GetSuperModuleNumber(vTowerIndex[cellIDLocal]), dEta, dPhi);
        }
      }
    } // end of bc loop
    for (const auto& trackIDs : vecCellToTracks) {
      cellMatchedTracks(trackIDs);
    }
  }
  PROCESS_SWITCH(EmcalCellTrackMatcher, processSortedWithSlices, "run analysis with sorted tracks utilizing preslices", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalCellTrackMatcher>(cfgc)};
}
