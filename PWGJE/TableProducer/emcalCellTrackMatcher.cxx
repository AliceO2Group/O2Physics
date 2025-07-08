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

#include "PWGJE/Core/EmcalMatchingUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DetectorsBase/GeometryManager.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/GeometryBase.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyGlobTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TracksDCA, o2::aod::TracksCov>;
using BcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredTracks = o2::soa::Filtered<MyGlobTracks>;
using McCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;
using FilteredCells = o2::soa::Filtered<aod::Calos>;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = bg::model::point<float, 2, bg::cs::cartesian>;
using box = bg::model::box<point>;
using value = std::pair<box, int>; // cell bounding box + cellID

struct EmcalCellTrackMatcher {

  // small neighbor cell struct
  struct Neighbor {
    int sm, iphi, ieta;
    bool valid;
    Neighbor() : sm(0), iphi(0), ieta(0), valid(true)
    {
    }
  };

  // small track data struct
  struct Track {
    float eta;
    float phi;
    int trackID;
  };

  Produces<o2::aod::EMCALCellTracks> cellMatchedTracks;

  // Options for the clusterization
  // 1 corresponds to EMCAL cells based on the Run2 definition.
  Configurable<int> selectedCellType{"selectedCellType", 1, "EMCAL Cell type"};
  Configurable<float> trackMinPt{"trackMinPt", 0.3f, "Minimum pT for tracks to perform track matching, to reduce computing time. Tracks below a certain pT will be loopers anyway."};
  Configurable<float> trackDCAz{"trackDCAz", 5.f, "Maximum DCAz of a track to to be considered primary and perform the track matching."};
  Configurable<float> trackDCAxy{"trackDCAxy", 5.f, "Maximum DCAxy of a track to to be considered primary and perform the track matching."};
  Configurable<float> maxMatchingDistance{"maxMatchingDistance", 0.01011162697, "Max matching distance track-cell. Default value is equal to half a cell diagonal."};
  Configurable<float> maxHalfDistance{"maxHalfDistance", 0.00715f, "Max distance in eta and phi between to cells. Default value is half a cell width, which is 0.0143."};
  Configurable<int> maxNumberTracks{"maxNumberTracks", 1000, "Number of tracks we expect per BC that point to EMCal. This is for optimization."};
  Configurable<bool> useAlignmentFromCCDB{"useAlignmentFromCCDB", true, "States if alignment objects should be used from CCDB"};
  Configurable<float> minAplitude{"minAplitude", 0.09f, "Minimum amplitude a cell needs to have to be matched with a track. Since we do not use cells below a certain threshold for clustering this can be used to speed up this task."};

  // Require EMCAL cells (CALO type 1)
  Partition<aod::Calos> emcalCells = aod::calo::caloType == selectedCellType.value && aod::calo::amplitude >= minAplitude.value;
  Filter emccellfilter = aod::calo::caloType == selectedCellType.value;

  // Filter for the tracks
  const float trackNotOnEMCal = -900;
  Filter trackFilter = (aod::track::pt >= trackMinPt) && (aod::track::trackEtaEmcal > trackNotOnEMCal) && (aod::track::trackPhiEmcal > trackNotOnEMCal) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (nabs(aod::track::dcaXY) < trackDCAxy) && (nabs(aod::track::dcaZ) < trackDCAz);

  // CDB service (for geometry)
  Service<o2::ccdb::BasicCCDBManager> mCcdbManager;

  // QA
  o2::framework::HistogramRegistry mHistManager{"EmcalCellTrackMatcherQAHistograms"};

  // EMCal geometry
  o2::emcal::Geometry* geometry;
  static constexpr int16_t NCellsInEMCal = 12288;
  static constexpr int NCellsTotal = 17664; // Number of cells in the EMCal
  static constexpr float Epsilon = 1.e-6;   // small value to ensure cells boxes are unique
  std::array<float, NCellsTotal> arrEta;
  std::array<float, NCellsTotal> arrPhi;
  std::array<float, NCellsTotal> arrLowEta;
  std::array<float, NCellsTotal> arrUpEta;
  std::array<float, NCellsTotal> arrLowPhi;
  std::array<float, NCellsTotal> arrUpPhi;

  // Preslices
  Preslice<MyGlobTracks> perCollision = o2::aod::track::collisionId;
  PresliceUnsorted<CollEventSels> collisionsPerFoundBC = aod::evsel::foundBCId;
  Preslice<aod::Calos> cellsPerFoundBC = aod::calo::bcId;
  Preslice<aod::SortedTracks> tracksPerFoundBC = aod::emcalcluster::bcId;

  Neighbor computeNeighbor(int16_t sm, int16_t iphi, int16_t ieta, int16_t dPhi, int16_t dEta, int16_t maxPhi, int16_t maxEta)
  {
    Neighbor neighbour;
    iphi += dPhi;
    ieta += dEta;

    if (iphi == maxPhi) {
      if (sm == 10 || sm == 11 || sm == 18 || sm == 19) { // o2-linter: disable=magic-number (these are SM that do not have another SM directly next to them in positive iphi direction)
        neighbour.valid = false;
      } else {
        iphi = 0;
        sm += 2;
      }
    } else if (iphi < 0) {
      if ((sm == 0 || sm == 1 || sm == 12 || sm == 13) || ((sm == 18 && ieta >= 32) || (sm == 19 && ieta <= 16))) { // o2-linter: disable=magic-number (these are SM that do not have another SM directly next to them in negative iphi direction. SM 18 and 19 only partially have a SM next to them depending on eta)
        neighbour.valid = false;
      } else {
        iphi = 23;
        sm -= 2;
      }
    }

    if (ieta == maxEta) {
      if ((sm % 2) || (sm == 12 || sm == 14 || sm == 16)) { // o2-linter: disable=magic-number (these are SM that do not have another SM directly next to them in positive ieta direction)
        neighbour.valid = false;
      } else {
        ieta = 0;
        sm += 1;
      }
    } else if (ieta < 0) {
      if ((sm % 2 == 0) || (sm == 13 || sm == 15 || sm == 17)) { // o2-linter: disable=magic-number (these are SM that do not have another SM directly next to them in negative ieta direction)
        neighbour.valid = false;
      } else {
        ieta = 47;
        sm -= 1;
      }
    }
    neighbour.ieta = ieta;
    neighbour.iphi = iphi;
    neighbour.sm = sm;
    return neighbour;
  };

  // Build the R-tree
  bgi::rtree<value, bgi::quadratic<16>> buildCellRTree(const std::vector<std::tuple<float, float, float, float>>& cellBounds)
  {
    std::vector<value> entries;
    int cellID = 0;
    for (const auto& [etaMin, etaMax, phiMin, phiMax] : cellBounds) {
      box b(point(etaMin, phiMin), point(etaMax, phiMax));
      entries.emplace_back(b, cellID++);
    }
    return bgi::rtree<value, bgi::quadratic<16>>(entries);
  }

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
    LOG(debug) << "Completed init!";

    // Setup QA hists.
    using O2HistType = o2::framework::HistType;
    o2::framework::AxisSpec etaAxis{1000, -0.8, 0.8, "#it{#eta}"},
      phiAxis{1000, 0, 2 * 3.14159, "#it{#varphi} (rad)"},
      dEtaAxis{100, -0.015, 0.015, "#Delta#it{#eta}"},
      dPhiAxis{100, -0.015, 0.015, "#Delta#it{#varphi} (rad)"},
      cellAxis{NCellsTotal, -0.5, NCellsTotal - 0.5, "cellID"},
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
    arrEta.fill(-999.f);
    arrPhi.fill(-999.f);
    for (int iCell = 0; iCell < NCellsTotal; ++iCell) {
      // get the local coordinates of the cell
      try {
        geometry->RelPosCellInSModule(iCell, dist).GetCoordinates(lxyzi[0], lxyzi[1], lxyzi[2]);
      } catch (o2::emcal::InvalidCellIDException& e) {
        LOG(error) << e.what();
        continue;
      }
      // Now get the global coordinate
      int iSM = geometry->GetSuperModuleNumber(iCell);
      geometry->GetGlobal(lxyzi, xyzi, iSM);
      auto [iRow, iCol] = geometry->GlobalRowColFromIndex(iCell);
      math_utils::Point3D<float> pos(xyzi[0], xyzi[1], xyzi[2]);
      arrEta[iCell] = pos.Eta();
      arrPhi[iCell] = RecoDecay::constrainAngle(pos.Phi());
      mHistManager.fill(HIST("hEMCalEtaPhi3DMap"), iCol, iRow, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusSM"), iSM, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusEta"), iCol, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
      mHistManager.fill(HIST("pEMCalRadiusPhi"), iRow, std::sqrt(pos.X() * pos.X() + pos.Y() * pos.Y()));
    }

    // second loop over the cells to get the cell neighbor values
    for (int iCell = 0; iCell < NCellsTotal; ++iCell) {
      auto [iSM, iMod, iIphi, iIeta] = geometry->GetCellIndex(iCell);
      auto [iphi, ieta] = geometry->GetCellPhiEtaIndexInSModule(iSM, iMod, iIphi, iIeta);

      int nRows = o2::emcal::EMCAL_ROWS; // rows -> phi
      int nCols = o2::emcal::EMCAL_COLS; // cols -> eta

      switch (geometry->GetSMType(iSM)) {
        case o2::emcal::EMCAL_STANDARD:
          // all good
          break;
        case o2::emcal::EMCAL_HALF:
          // does not exist
          break;
        case o2::emcal::DCAL_STANDARD:
          nCols = 2 * nCols / 3;
          break;
        case o2::emcal::EMCAL_THIRD:
        case o2::emcal::DCAL_EXT:
          nRows = nRows / 3;
          break;
        default:
          break;
      }

      Neighbor nPhiUp = computeNeighbor(iSM, iphi, ieta, +1, 0, nRows, nCols);
      Neighbor nPhiDown = computeNeighbor(iSM, iphi, ieta, -1, 0, nRows, nCols);
      Neighbor nEtaUp = computeNeighbor(iSM, iphi, ieta, 0, +1, nRows, nCols);
      Neighbor nEtaDown = computeNeighbor(iSM, iphi, ieta, 0, -1, nRows, nCols);

      auto computeBoundary = [&](float center, const Neighbor& n, bool isEta, bool positiveDir) -> float {
        if (!n.valid)
          return center + (positiveDir ? +maxHalfDistance.value : -maxHalfDistance.value);

        int neighId = geometry->GetAbsCellIdFromCellIndexes(n.sm, n.iphi, n.ieta);
        float neighVal = isEta ? arrEta[neighId] : arrPhi[neighId];
        float delta = std::abs(center - neighVal) / 2.f - Epsilon;
        if (delta >= maxHalfDistance.value) {
          delta = maxHalfDistance.value;
        }

        return center + (positiveDir ? +delta : -delta);
      };
      arrUpPhi[iCell] = computeBoundary(arrPhi[iCell], nPhiUp, false, true);
      arrLowPhi[iCell] = computeBoundary(arrPhi[iCell], nPhiDown, false, false);
      arrUpEta[iCell] = computeBoundary(arrEta[iCell], nEtaUp, true, true);
      arrLowEta[iCell] = computeBoundary(arrEta[iCell], nEtaDown, true, false);

    } // for (int iCell = 0; iCell < NCellsTotal; ++iCell) {
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
              if (cellNumber < 0 || cellNumber >= NCellsTotal) {
                continue;
              }
              vCellEta.emplace_back(arrEta[cellNumber]);
              vCellPhi.emplace_back(arrPhi[cellNumber]);
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
                if (NCellsInEMCal > vTowerIndex[cellIDLocal]) {
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
    cellMatchedTracks.reserve(vecCellToTracks.size());
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

    // if we have no cells, just skip this DF
    if (!cells.size()) {
      return;
    }
    // if we have no tracks fill the table with empty spans
    if (!tracks.size()) {
      cellMatchedTracks.reserve(cells.size());
      for (int64_t iCell = 0; iCell < cells.size(); ++iCell) {
        cellMatchedTracks(std::span<int>());
      }
      return;
    }

    // outer vector is size as number of ALL calo cells in current DF, inside vector will store trackIDs
    std::vector<std::vector<int>> vecCellToTracks(cells.size(), std::vector<int>());
    // std::vector<int> vecTrackIDs;
    // vecTrackIDs.reserve(cells.size());
    // std::vector<int> vecNTrackIDs(cells.size(), 0);

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
    vCellPhi.reserve(NCellsTotal);
    vCellEta.reserve(NCellsTotal);
    vCellIndex.reserve(NCellsTotal);
    vTowerIndex.reserve(NCellsTotal);
    vSMIndex.reserve(NCellsTotal);
    while (cell != cellEnd && track != trackEnd) {
      auto cellBC = cell.bcId();
      auto trackBC = track.bcId();

      if (cellBC == trackBC) {
        auto currentBCID = cellBC;

        while (cell != cellEnd && cellBC == currentBCID) {
          int16_t cellNumber = cell.cellNumber();
          if (cellNumber < 0 || cellNumber >= NCellsTotal) {
            LOG(info) << "cell number " << cellNumber << " not within EMCal!";
            continue;
          }
          vCellEta.emplace_back(arrEta[cellNumber]);
          vCellPhi.emplace_back(arrPhi[cellNumber]);
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
            // vecTrackIDs[cellID]
            // vecNTrackIDs[cellID]
            mHistManager.fill(HIST("hvTrackEtaPhiEMCal"), trackEta, trackPhi);
            float dEta = trackEta - vCellEta[cellIDLocal];
            float dPhi = trackPhi - vCellPhi[cellIDLocal];
            if (NCellsInEMCal > vTowerIndex[cellIDLocal]) {
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

    // vecTrackIDs
    // vecNTrackIDs

    cellMatchedTracks.reserve(vecCellToTracks.size());
    for (const auto& trackIDs : vecCellToTracks) {
      cellMatchedTracks(trackIDs);
    }
    LOG(info) << "number of entries in new matched tracks table " << cellMatchedTracks.lastIndex() << "\t number of inital cell table entries " << cells.size();
  }
  PROCESS_SWITCH(EmcalCellTrackMatcher, processSorted, "run analysis with tracks sorted according to BCId", false);

  void processSortedRTree(o2::aod::SortedTracks const& tracks, aod::Calos const& cells)
  {
    LOG(debug) << "Starting process sorted tracks.";

    // if we have no cells, just skip this DF
    if (cells.size() == 0) {
      return;
    }

    auto emptyVec = std::vector<int>({});
    auto emptySpan = std::span(emptyVec.begin(), 0);
    // if we have no tracks fill the table with empty spans
    cellMatchedTracks.reserve(cells.size());
    if (tracks.size() == 0) {
      for (int64_t iCell = 0; iCell < cells.size(); ++iCell) {
        cellMatchedTracks(emptySpan);
      }
      return;
    } else {
      for (int64_t iCell = 0; iCell < cells.size(); ++iCell) {
        cellMatchedTracks(emptySpan);
      }
    }

    // outer vector is size as number of ALL calo cells in current DF, inside vector will store trackIDs
    // std::vector<std::vector<int>> vecCellToTracks(cells.size(), std::vector<int>());
    // std::vector<int> vecTrackIDs;
    // vecTrackIDs.reserve(cells.size());
    // std::vector<int> vecNTrackIDs(cells.size(), 0);

    // std::array<size_t, NCellsTotal> matchedTrackCounts;
    // std::vector<size_t> matchedTrackOffsets;
    // matchedTrackOffsets.reserve(NCellsTotal + 1);
    // std::vector<int> matchedTrackCountsOrdered;
    // matchedTrackCountsOrdered.reserve(NCellsTotal);

    // auto track = tracks.begin();
    // auto cell = emcalCells.begin();

    // const auto cellEnd = emcalCells.end();
    // const auto trackEnd = tracks.end();

    // std::vector<int> matchIndexCell;
    // matchIndexCell.reserve(maxNumberTracks); // reserve enough space for better performance
    // std::vector<int> matchIndexTower;
    // matchIndexTower.reserve(maxNumberTracks); // reserve enough space for better performance

    // // For the matching of a track to a cell we will use a KDTree approach from emcmatchingutilities
    // // where cells are matched to tracks. For this we need vectors of cell and track eta and phi.
    // std::vector<float> vTrackPhi;
    // std::vector<float> vTrackEta;
    // std::vector<int> vTrackIndex;
    // std::vector<float> vCellPhi;
    // std::vector<float> vCellEta;
    // std::vector<int> vCellIndex;
    // std::vector<int16_t> vTowerIndex;
    // std::vector<int> vSMIndex;
    // std::vector<std::tuple<float, float, float, float>> cellBounds;
    // vTrackPhi.reserve(maxNumberTracks);
    // vTrackEta.reserve(maxNumberTracks);
    // vTrackIndex.reserve(maxNumberTracks);
    // vCellPhi.reserve(NCellsTotal);
    // vCellEta.reserve(NCellsTotal);
    // vCellIndex.reserve(NCellsTotal);
    // vTowerIndex.reserve(NCellsTotal);
    // vSMIndex.reserve(NCellsTotal);
    // cellBounds.reserve(NCellsTotal);
    // while (cell != cellEnd && track != trackEnd) {
    //   auto cellBC = cell.bcId();
    //   auto trackBC = track.bcId();

    //   if (cellBC == trackBC) {
    //     auto currentBCID = cellBC;

    //     while (cell != cellEnd && cellBC == currentBCID) {
    //       int16_t cellNumber = cell.cellNumber();
    //       if (cellNumber < 0 || cellNumber >= NCellsTotal) {
    //         LOG(info) << "cell number " << cellNumber << " not within EMCal!";
    //         continue;
    //       }
    //       vCellEta.emplace_back(arrEta[cellNumber]);
    //       vCellPhi.emplace_back(arrPhi[cellNumber]);
    //       vCellIndex.emplace_back(cell.globalIndex());
    //       vTowerIndex.emplace_back(cellNumber);
    //       vSMIndex.emplace_back(geometry->GetSuperModuleNumber(cellNumber));
    //       // cell bounds (etaMin, etaMax, phiMin, phiMax)
    //       cellBounds.emplace_back(std::tuple<float, float, float, float>(arrLowEta[cellNumber], arrUpEta[cellNumber], arrLowPhi[cellNumber], arrUpPhi[cellNumber]));
    //       if (++cell != cellEnd) {
    //         cellBC = cell.bcId();
    //       }
    //     }

    //     while (track != trackEnd && trackBC == currentBCID) {
    //       vTrackPhi.emplace_back(track.phi());
    //       vTrackEta.emplace_back(track.eta());
    //       vTrackIndex.emplace_back(track.trackId());
    //       if (++track != trackEnd) {
    //         trackBC = track.bcId();
    //       }
    //     }
    //     // build the current RTree
    //     auto cellRTree = buildCellRTree(cellBounds);

    //     // loop over tracks to find matching cell
    //     for (size_t iTrack = 0; iTrack < vTrackIndex.size(); ++iTrack) {
    //       point p(vTrackEta[iTrack], vTrackPhi[iTrack]);
    //       std::vector<value> result;
    //       cellRTree.query(bgi::contains(p), std::back_inserter(result));

    //       if (!result.empty()) {
    //         // Choose best match if multiple
    //         auto best = std::min_element(result.begin(), result.end(),
    //                                      [&](const value& a, const value& b) {
    //                                        auto center = [](const value& v) {
    //                                          float etaC = 0.5f * (v.first.min_corner().get<0>() + v.first.max_corner().get<0>());
    //                                          float phiC = 0.5f * (v.first.min_corner().get<1>() + v.first.max_corner().get<1>());
    //                                          return std::make_pair(etaC, phiC);
    //                                        };
    //                                        auto [etaA, phiA] = center(a);
    //                                        auto [etaB, phiB] = center(b);
    //                                        float dA = std::hypot(vTrackEta[iTrack] - etaA, vTrackPhi[iTrack] - phiA);
    //                                        float dB = std::hypot(vTrackEta[iTrack] - etaB, vTrackPhi[iTrack] - phiB);
    //                                        return dA < dB;
    //                                      });

    //         int localCellIdx = best->second;        // Index in vCellIndex
    //         int cellID = vTowerIndex[localCellIdx]; // towerID of the matched cell
    //         LOG(info) << "localCellIdx = " << localCellIdx << "\t cellID = " << cellID;
    //         matchIndexCell.emplace_back(localCellIdx);
    //         matchIndexTower.emplace_back(cellID);

    //         // Fill histograms
    //         mHistManager.fill(HIST("hvTrackEtaPhiEMCal"), vTrackEta[iTrack], vTrackPhi[iTrack]);
    //         float dEta = vTrackEta[iTrack] - vCellEta[localCellIdx];
    //         float dPhi = vTrackPhi[iTrack] - vCellPhi[localCellIdx];

    //         if (NCellsInEMCal > vTowerIndex[localCellIdx]) {
    //           mHistManager.fill(HIST("hTrackCellDiffEMCal"), dEta, dPhi);
    //         } else {
    //           mHistManager.fill(HIST("hTrackCellDiffDCal"), dEta, dPhi);
    //         }

    //         auto [iRow, iCol] = geometry->GlobalRowColFromIndex(vTowerIndex[localCellIdx]);
    //         mHistManager.fill(HIST("hNTracksPerCell"), iCol, iRow);
    //         fillSupermoduleHistograms(vSMIndex[localCellIdx], dEta, dPhi);

    //       } else {
    //         matchIndexCell.emplace_back(-1);
    //         matchIndexTower.emplace_back(-1);
    //       } // if (!result.empty())
    //     } // end of loop over tracks

    //     // count how many tracks a single cell is matched to
    //     for (size_t iTrack = 0; iTrack < matchIndexTower.size(); ++iTrack) {
    //       int towerID = matchIndexTower[iTrack];
    //       if (towerID != -1) {
    //         matchedTrackCounts[towerID]++;
    //       }
    //     }
    //     matchedTrackOffsets.assign(vTowerIndex.size() + 1, 0);
    //     matchedTrackCountsOrdered.assign(vTowerIndex.size(), 0);
    //     // make adjustments for the offset for the flat vector later
    //     for (size_t iCell = 0; iCell < vTowerIndex.size(); ++iCell) {
    //       matchedTrackOffsets[iCell + 1] = matchedTrackOffsets[iCell] + matchedTrackCounts[vTowerIndex[iCell]];
    //       matchedTrackCountsOrdered[iCell] = matchedTrackCounts[vTowerIndex[iCell]];
    //     }

    //     // making a flat vector order by the appearance of the cells
    //     std::vector<int> matchedTrackIDsFlat(matchedTrackOffsets.back());
    //     matchedTrackIDsFlat.assign(matchedTrackOffsets.back(), -1);
    //     std::vector<size_t> cellFillPos = matchedTrackOffsets; // working positions
    //     if (!matchedTrackIDsFlat.empty()) {
    //       for (size_t iTrack = 0; iTrack < matchIndexCell.size(); ++iTrack) {
    //         int localCellID = matchIndexCell[iTrack];
    //         if (localCellID != -1) {
    //           size_t pos = cellFillPos[localCellID]++;
    //           matchedTrackIDsFlat[pos] = vTrackIndex[iTrack];
    //         }
    //       }
    //     }
    //     auto it = matchedTrackIDsFlat.begin();
    //     LOG(info) << "matchedTrackCountsOrdered.size() = " << matchedTrackCountsOrdered.size();
    //     for (auto& size : matchedTrackCountsOrdered) {
    //       auto sp = std::span(it, size);
    //       // fill sp
    //       cellMatchedTracks(sp);
    //       it += size;
    //     }

    //     matchedTrackCounts.fill(0);
    //     matchedTrackCountsOrdered.clear();
    //     cellBounds.clear();
    //     matchedTrackOffsets.clear();
    //     vTrackPhi.clear();
    //     vTrackEta.clear();
    //     vTrackIndex.clear();
    //     vCellPhi.clear();
    //     vCellEta.clear();
    //     vCellIndex.clear();
    //     vTowerIndex.clear();
    //     vSMIndex.clear();
    //     matchIndexCell.clear();
    //     matchIndexTower.clear();
    //   } else if (cellBC < trackBC) {
    //     LOG(info) << "filling empty spans";
    //     while (cell != cellEnd && cellBC < trackBC) {
    //       if (++cell != cellEnd) {
    //         cellMatchedTracks(std::span<const int>());
    //         cellBC = cell.bcId();
    //       }
    //     }
    //   } else {
    //     while (track != trackEnd && trackBC < cellBC) {
    //       if (++track != trackEnd) {
    //         trackBC = track.bcId();
    //       }
    //     }
    //   }
    //   if (cell == cellEnd || track == trackEnd) {
    //     break;
    //   }
    // } // while (cell != cellEnd && track != trackEnd)

    // LOG(info) << "number of entries in new matched tracks table " << cellMatchedTracks.lastIndex() + 1 << "\t number of inital cell table entries " << cells.size();
  }
  PROCESS_SWITCH(EmcalCellTrackMatcher, processSortedRTree, "run analysis with tracks sorted according to BCId", false);

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
        vCellEta.emplace_back(arrEta[cellNumber]);
        vCellPhi.emplace_back(arrPhi[cellNumber]);
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
          if (NCellsInEMCal > vTowerIndex[cellIDLocal]) {
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
    cellMatchedTracks.reserve(vecCellToTracks.size());
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
