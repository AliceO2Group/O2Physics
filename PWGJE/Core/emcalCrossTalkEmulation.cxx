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

/// \file emcalCrossTalkEmulation.cxx
/// \brief emulation of emcal cross talk for simulations
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe-University

#include "emcalCrossTalkEmulation.h"

#include <DataFormatsEMCAL/Cell.h>
#include <DataFormatsEMCAL/CellLabel.h>
#include <DataFormatsEMCAL/Constants.h>
#include <EMCALBase/Geometry.h>
#include <Framework/Array2D.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <algorithm> // std::find_if
#include <array>
#include <cstddef> // size_t
#include <cstdlib> // std::abs
#include <iomanip> // setw
#include <ios>     // left and right
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
// #include "Framework/OutputObjHeader.h"

// #include "Common/CCDB/EventSelectionParams.h"
#include <TH1.h>
#include <TRandom3.h>

#include <fairlogger/Logger.h>

using namespace o2;
using namespace o2::emccrosstalk;
using namespace o2::framework;

template <typename T, std::size_t N>
auto printArray(std::array<T, N> const& arr)
{
  std::stringstream ss;
  ss << "\n[SM0: " << arr[0];
  for (auto i = 1u; i < N; ++i) {
    ss << ", SM" << i << ": " << arr[i];
  }
  ss << "]";
  return ss.str();
}

template <typename T>
auto printMatrix(Array2D<T> const& m)
{
  std::stringstream ss;
  // Print column headers
  ss << std::endl
     << std::setw(6) << " " << std::setw(10) << "value1"
     << std::setw(10) << "value2"
     << std::setw(10) << "value3"
     << std::setw(10) << "value4" << std::endl;

  // Print rows with SM labels
  for (size_t i = 0; i < m.rows; ++i) {
    ss << "SM" << std::left << std::setw(3) << i; // e.g., SM0, SM1...
    for (size_t j = 0; j < m.cols; ++j) {
      ss << std::right << std::setw(10) << m(i, j);
    }
    ss << std::endl;
  }

  return ss.str();
}

void init2DElement(Array2D<float>& matrix, const Array2D<float>& config, const char* name)
{
  int rows = config.rows;
  int cols = config.cols;

  if (rows == 0 && cols == 0) {
    LOG(info) << name << " has size 0 x 0, so it is disabled!";
  } else if (cols != NNeighbourCases || (rows != 1 && rows != NSM)) {
    LOG(error) << name << " must have 4 columns and either 1 or 20 rows!";
  } else {
    for (int sm = 0; sm < NSM; ++sm) {
      const int row = (rows == 1) ? 0 : sm;

      for (int i = 0; i < cols; ++i) {
        matrix[sm][i] = config(row, i);
      }
    }
  }
}

void init1DElement(std::array<float, NSM>& arr, const std::vector<float>& config, const char* name)
{
  size_t confSize = config.size();
  if (confSize == 0) {
    LOG(info) << name << " has size 0, so it is disabled!";
  } else if (config.size() != 1 && confSize != NSM) {
    LOG(error) << name << " must have either size 1 or 20!";
  } else {
    for (int sm = 0; sm < NSM; ++sm) {
      const int row = (confSize == 1) ? 0 : sm;
      arr[sm] = config[row];
    }
  }
}

o2::framework::AxisSpec axisEnergy = {7000, 0.f, 70.f, "#it{E}_{cell} (GeV)"};
// For each of the following configurables we will use:
// empty vector == disabled
// vector of size 4 == same for all SM
// vector of vectors with size nSM * 4 == each SM has its own setting
// the 4 values (0-3) correspond to in relative [row,col]: 0: [+-1,0], 1: [+-1,+or-1], 2:  [0,+or-1], 3: [+-2, 0 AND +or-1]
// +---+---+-----+---+---+--+--+--+
// | 3 | 0 | Hit | 0 | 3 |  |  |  |
// +---+---+-----+---+---+--+--+--+
// | 3 | 1 |  2  | 1 | 3 |  |  |  |
// +---+---+-----+---+---+--+--+--+

void EMCCrossTalk::initObjects(const EmcCrossTalkConf& config)
{
  const int run3RunNumber = 223409;
  mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(run3RunNumber);
  if (!mGeometry) {
    LOG(error) << "Failure accessing mGeometry";
  }

  // first set the simple run time variables
  mTCardCorrClusEnerConserv = config.conserveEnergy.value;
  mRandomizeTCard = config.randomizeTCardInducedEnergy.value;
  mTCardCorrMinAmp = config.inducedTCardMinimumCellEnergy.value;
  mTCardCorrMinInduced = config.inducedTCardMinimum.value;
  mTCardCorrMaxInducedELeak = config.inducedTCardMaximumELeak.value;
  mTCardCorrMaxInduced = config.inducedTCardMaximum.value;

  // 2nd define the NSM x NNeighbourCases matrices
  mTCardCorrInduceEner = Array2D<float>(std::vector<float>(NSM * 4, 0.f), NSM, 4);
  mTCardCorrInduceEnerFrac = Array2D<float>(std::vector<float>(NSM * 4, 0.f), NSM, 4);
  mTCardCorrInduceEnerFracP1 = Array2D<float>(std::vector<float>(NSM * 4, 0.f), NSM, 4);
  mTCardCorrInduceEnerFracWidth = Array2D<float>(std::vector<float>(NSM * 4, 0.f), NSM, 4);

  // now properly init the NSM x NNeighbourCases matrices
  // ------------------------------------------------------------------------
  // mTCardCorrInduceEner
  init2DElement(mTCardCorrInduceEner, config.inducedEnergyLossConstant.value, "inducedEnergyLossConstant");

  // mTCardCorrInduceEnerFrac
  init2DElement(mTCardCorrInduceEnerFrac, config.inducedEnergyLossFraction.value, "inducedEnergyLossFraction");

  // mTCardCorrInduceEnerFracP1
  init2DElement(mTCardCorrInduceEnerFracP1, config.inducedEnergyLossFractionP1.value, "inducedEnergyLossFractionP1");

  // mTCardCorrInduceEnerFracWidth
  init2DElement(mTCardCorrInduceEnerFracWidth, config.inducedEnergyLossFractionWidth.value, "inducedEnergyLossFractionWidth");
  // ------------------------------------------------------------------------

  init1DElement(mTCardCorrInduceEnerFracMax, config.inducedEnergyLossMaximumFraction.value, "inducedEnergyLossMaximumFraction");
  init1DElement(mTCardCorrInduceEnerFracMin, config.inducedEnergyLossMinimumFraction.value, "inducedEnergyLossMinimumFraction");
  init1DElement(mTCardCorrInduceEnerFracMinCentralEta, config.inducedEnergyLossMinimumFractionCentralEta.value, "inducedEnergyLossMinimumFractionCentralEta");
  init1DElement(mTCardCorrInduceEnerProb, config.inducedEnergyLossProbability.value, "inducedEnergyLossProbability");

  resetArrays();

  // Print the full matrices and vectors that will be used:
  if (config.printConfiguration.value) {
    LOGF(info, "inducedEnergyLossConstant: %s", printMatrix((mTCardCorrInduceEner)));
    LOGF(info, "inducedEnergyLossFraction: %s", printMatrix((mTCardCorrInduceEnerFrac)));
    LOGF(info, "inducedEnergyLossFractionP1: %s", printMatrix((mTCardCorrInduceEnerFracP1)));
    LOGF(info, "inducedEnergyLossFractionWidth: %s", printMatrix((mTCardCorrInduceEnerFracWidth)));
    LOGF(info, "inducedEnergyLossMaximumFraction: %s", printArray(mTCardCorrInduceEnerFracMax).c_str());
    LOGF(info, "inducedEnergyLossMinimumFraction: %s", printArray(mTCardCorrInduceEnerFracMin).c_str());
    LOGF(info, "inducedEnergyLossMinimumFractionCentralEta: %s", printArray(mTCardCorrInduceEnerFracMinCentralEta).c_str());
    LOGF(info, "inducedEnergyLossProbability: %s", printArray(mTCardCorrInduceEnerProb).c_str());
  }
}

void EMCCrossTalk::resetArrays()
{
  for (size_t j = 0; j < NCells; j++) {
    mTCardCorrCellsEner[j] = 0.;
    mTCardCorrCellsNew[j] = false;
  }

  mCellsTmp.clear();
}

void EMCCrossTalk::setCells(std::vector<o2::emcal::Cell>& cells, std::vector<o2::emcal::CellLabel>& cellLabels)
{
  mCells = &cells;
  mCellsTmp = cells; // a copy since we will need one vector with the changed energies and one with the original ones
  mCellLabels = &cellLabels;
}

void EMCCrossTalk::calculateInducedEnergyInTCardCell(int absId, int absIdRef, int iSM, float ampRef, int cellCase)
{
  // Check that the cell exists
  if (absId < 0) {
    return;
  }

  // Get the fraction
  float frac = mTCardCorrInduceEnerFrac[iSM][cellCase] + ampRef * mTCardCorrInduceEnerFracP1[iSM][cellCase];

  // Use an absolute minimum and maximum fraction if calculated one is out of range
  if (frac < mTCardCorrInduceEnerFracMin[iSM]) {
    frac = mTCardCorrInduceEnerFracMin[iSM];
  } else if (frac > mTCardCorrInduceEnerFracMax[iSM]) {
    frac = mTCardCorrInduceEnerFracMax[iSM];
  }

  // If active, use different absolute minimum fraction for central eta, exclude DCal 2/3 SM
  if (mTCardCorrInduceEnerFracMinCentralEta[iSM] > 0 && (iSM < FirstDCal23SM || iSM > LastDCal23SM)) {
    // Odd  SM
    int ietaMin = 32;
    int ietaMax = 47;

    // Even SM
    if (iSM % 2) {
      ietaMin = 0;
      ietaMax = 15;
    }

    // First get the SM, col-row of this tower
    // int imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
    auto [iSM, iMod, iIphi, iIeta] = mGeometry->GetCellIndex(absId);
    auto [iphi, ieta] = mGeometry->GetCellPhiEtaIndexInSModule(iSM, iMod, iIphi, iIeta);

    if (ieta >= ietaMin && ieta <= ietaMax) {
      if (frac < mTCardCorrInduceEnerFracMinCentralEta[iSM])
        frac = mTCardCorrInduceEnerFracMinCentralEta[iSM];
    }
  } // central eta

  LOGF(debug, "\t fraction %2.3f", frac);

  // Randomize the induced fraction, if requested
  if (mRandomizeTCard) {
    frac = mRandom.Gaus(frac, mTCardCorrInduceEnerFracWidth[iSM][cellCase]);
    LOGF(debug, "\t randomized fraction %2.3f", frac);
  }

  // If fraction too small or negative, do nothing else
  if (frac < Epsilon) {
    return;
  }

  // Calculate induced energy
  float inducedE = mTCardCorrInduceEner[iSM][cellCase] + ampRef * frac;

  // Check if we induce too much energy, in such case use a constant value
  if (mTCardCorrMaxInduced < inducedE)
    inducedE = mTCardCorrMaxInduced;

  LOGF(debug, "\t induced E %2.3f", inducedE);

  // Try to find the cell that will get energy induced
  float amp = 0.f;
  auto itCell = std::find_if((*mCells).begin(), (*mCells).end(), [absId](const o2::emcal::Cell& cell) {
    return cell.getTower() == absId;
  });

  if (itCell != (*mCells).end()) {
    // We found a cell, so let's get the amplitude of that cell
    amp = itCell->getAmplitude();
  } else {
    amp = 0.f; // this is a new cell, so the base amp is 0.f
  }

  // Check that the induced+amp is large enough to avoid extra linearity effects
  // typically of the order of the clusterization cell energy cut
  // if inducedTCardMaximumELeak was set to a positive value, then induce the energy as long as its smaller than that value
  if ((amp + inducedE) > mTCardCorrMinInduced || inducedE < mTCardCorrMaxInducedELeak) {
    mTCardCorrCellsEner[absId] += inducedE;

    // If original energy of cell was null, create new one
    if (amp <= Epsilon) {
      mTCardCorrCellsNew[absId] = true;
    }
  } else {
    return;
  }

  LOGF(debug, "Cell %d is with amplitude %2.3f GeV is inducing %1.3f GeV energy to cell %d which already has %2.3f GeV energy with fraction %1.5f", absIdRef, ampRef, inducedE, absId, amp, frac);

  // Subtract the added energy to main cell, if energy conservation is requested
  if (mTCardCorrClusEnerConserv) {
    mTCardCorrCellsEner[absIdRef] -= inducedE;
  }
}

void EMCCrossTalk::makeCellTCardCorrelation()
{
  int id = -1;
  float amp = -1;

  // Loop on all cells with signal
  for (const auto& cell : (*mCells)) {
    id = cell.getTower();
    amp = cell.getAmplitude();

    if (amp <= mTCardCorrMinAmp) {
      continue;
    }

    // First get the SM, col-row of this tower
    auto [iSM, iMod, iIphi, iIeta] = mGeometry->GetCellIndex(id);
    auto [iphi, ieta] = mGeometry->GetCellPhiEtaIndexInSModule(iSM, iMod, iIphi, iIeta);

    // Determine randomly if we want to create a correlation for this cell,
    // depending the SM number of the cell
    if (mTCardCorrInduceEnerProb[iSM] < 1) {
      if (mRandom.Uniform(0, 1) > mTCardCorrInduceEnerProb[iSM]) {
        continue;
      }
    }

    LOGF(debug, "Reference cell absId %d, iEta %d, iPhi %d, amp %2.3f", id, ieta, iphi, amp);

    // Get the absId of the cells in the cross and same T-Card
    int absIDup = -1;
    int absIDdo = -1;
    int absIDlr = -1;
    int absIDuplr = -1;
    int absIDdolr = -1;

    int absIDup2 = -1;
    int absIDup2lr = -1;
    int absIDdo2 = -1;
    int absIDdo2lr = -1;

    // Only 2 columns in the T-Card, +1 for even and -1 for odd with respect reference cell
    // Sine we only have full T-Cards, we do not need to make any edge case checks
    // There is always either a column (eta direction) below or above
    int colShift = +1;
    if (ieta % 2) {
      colShift = -1;
    }

    absIDlr = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi, ieta + colShift);

    // Check if up / down cells from reference cell are not out of SM
    // First check if there is space one above
    if (iphi < emcal::EMCAL_ROWS - 1) {
      absIDup = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi + 1, ieta);
      absIDuplr = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi + 1, ieta + colShift);
    }

    // 2nd check if there is space one below
    if (iphi > 0) {
      absIDdo = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi - 1, ieta);
      absIDdolr = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi - 1, ieta + colShift);
    }

    // 3rd check if there is space two above
    if (iphi < emcal::EMCAL_ROWS - 2) {
      absIDup2 = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi + 2, ieta);
      absIDup2lr = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi + 2, ieta + colShift);
    }

    // 4th check if there is space two below
    if (iphi > 1) {
      absIDdo2 = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi - 2, ieta);
      absIDdo2lr = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphi - 2, ieta + colShift);
    }

    // Check if those cells are in the same T-Card
    int tCard = iphi / 8;
    if (tCard != (iphi + 1) / 8) {
      absIDup = -1;
      absIDuplr = -1;
    }
    if (tCard != (iphi - 1) / 8) {
      absIDdo = -1;
      absIDdolr = -1;
    }
    if (tCard != (iphi + 2) / 8) {
      absIDup2 = -1;
      absIDup2lr = -1;
    }
    if (tCard != (iphi - 2) / 8) {
      absIDdo2 = -1;
      absIDdo2lr = -1;
    }

    // Calculate induced energy to T-Card cells
    // first check if for the given cell case we actually do induce some energy
    if (((std::abs(mTCardCorrInduceEner[iSM][0]) > Epsilon) || (std::abs(mTCardCorrInduceEnerFrac[iSM][0]) > Epsilon)) && (std::abs(mTCardCorrInduceEnerFracP1[iSM][0]) > Epsilon) && (std::abs(mTCardCorrInduceEnerFracWidth[iSM][0]) > Epsilon)) {
      if (absIDup >= 0) {
        LOGF(debug, "cell up %d:", absIDup);
        calculateInducedEnergyInTCardCell(absIDup, id, iSM, amp, 0);
      }
      if (absIDdo >= 0) {
        LOGF(debug, "cell down %d:", absIDdo);
        calculateInducedEnergyInTCardCell(absIDdo, id, iSM, amp, 0);
      }
    }
    if (((std::abs(mTCardCorrInduceEner[iSM][1]) > Epsilon) || (std::abs(mTCardCorrInduceEnerFrac[iSM][1]) > Epsilon)) && (std::abs(mTCardCorrInduceEnerFracP1[iSM][1]) > Epsilon) && (std::abs(mTCardCorrInduceEnerFracWidth[iSM][1]) > Epsilon)) {
      if (absIDuplr >= 0) {
        LOGF(debug, "cell up left-right %d:", absIDuplr);
        calculateInducedEnergyInTCardCell(absIDuplr, id, iSM, amp, 1);
      }
      if (absIDdolr >= 0) {
        LOGF(debug, "cell down left-right %d:", absIDdolr);
        calculateInducedEnergyInTCardCell(absIDdolr, id, iSM, amp, 1);
      }
    }
    if (((std::abs(mTCardCorrInduceEner[iSM][2]) > Epsilon) || (std::abs(mTCardCorrInduceEnerFrac[iSM][2]) > Epsilon)) && (std::abs(mTCardCorrInduceEnerFracP1[iSM][2]) > Epsilon) && (std::abs(mTCardCorrInduceEnerFracWidth[iSM][2]) > Epsilon)) {
      if (absIDlr >= 0) {
        LOGF(debug, "cell left-right %d:", absIDlr);
        calculateInducedEnergyInTCardCell(absIDlr, id, iSM, amp, 2);
      }
    }
    if (((std::abs(mTCardCorrInduceEner[iSM][3]) > Epsilon) || (std::abs(mTCardCorrInduceEnerFrac[iSM][3]) > Epsilon)) && (std::abs(mTCardCorrInduceEnerFracP1[iSM][3]) > Epsilon) && (std::abs(mTCardCorrInduceEnerFracWidth[iSM][3]) > Epsilon)) {
      if (absIDup2 >= 0) {
        LOGF(debug, "cell up 2nd row %d:", absIDup2);
        calculateInducedEnergyInTCardCell(absIDup2, id, iSM, amp, 3);
      }
      if (absIDdo2 >= 0) {
        LOGF(debug, "cell down 2nd row %d:", absIDdo2);
        calculateInducedEnergyInTCardCell(absIDdo2, id, iSM, amp, 3);
      }
      if (absIDup2lr >= 0) {
        LOGF(debug, "cell up left-right 2nd row %d:", absIDup2lr);
        calculateInducedEnergyInTCardCell(absIDup2lr, id, iSM, amp, 3);
      }
      if (absIDdo2lr >= 0) {
        LOGF(debug, "cell down left-right 2nd row %d:", absIDdo2lr);
        calculateInducedEnergyInTCardCell(absIDdo2lr, id, iSM, amp, 3);
      }
    }
  } // cell loop
}

void EMCCrossTalk::addInducedEnergiesToExistingCells()
{
  // Add the induced energy to the cells and copy them into a new temporal container
  // used in AddInducedEnergiesToNewCells() to refill the default cells list fCaloCells
  // Create the data member only once. Done here, not sure where to do this properly in the framework.

  for (auto& cell : mCellsTmp) { // o2-linter: disable=const-ref-in-for-loop (we are changing a value here)
    float amp = cell.getAmplitude() + mTCardCorrCellsEner[cell.getTower()];
    // Set new amplitude in new temporal container
    cell.setAmplitude(amp);
  }
}

void EMCCrossTalk::addInducedEnergiesToNewCells()
{
  // count how many new cells
  size_t nCells = (*mCells).size();
  int nCellsNew = 0;
  for (size_t j = 0; j < NCells; j++) {
    // Newly created?
    if (!mTCardCorrCellsNew[j]) {
      continue;
    }
    // Accept only if at least 10 MeV
    if (mTCardCorrCellsEner[j] < MinCellEnergy) {
      continue;
    }
    nCellsNew++;
  }

  // reserve more space for new cell entries in original cells and celllabels
  (*mCells).reserve(nCells + nCellsNew);
  (*mCellLabels).reserve(nCells + nCellsNew);

  // change the amplitude of the original cells using
  for (size_t iCell = 0; iCell < mCellsTmp.size(); ++iCell) {
    (*mCells)[iCell].setAmplitude(mCellsTmp[iCell].getAmplitude());
  }

  // Add the new cells
  int absId = -1;
  float amp = -1;
  float time = 0;
  std::vector<int32_t> mclabel;

  for (size_t j = 0; j < NCells; j++) {
    // Newly created?
    if (!mTCardCorrCellsNew[j]) {
      continue;
    }

    // Accept only if at least 10 MeV
    if (mTCardCorrCellsEner[j] < MinCellEnergy) {
      continue;
    }

    // Add new cell
    absId = j;
    amp = mTCardCorrCellsEner[j];
    time = 615.e-9f;
    mclabel = {-1};

    // Assign as MC label the label of the neighboring cell with highest energy
    // within the same T-Card. Follow same approach for time.
    // Simplest assumption, not fully correct.
    // Still assign 0 as fraction of energy.

    // First get the iphi and ieta of this tower
    auto [iSM, iMod, iIphi, iIeta] = mGeometry->GetCellIndex(absId);
    auto [iphi, ieta] = mGeometry->GetCellPhiEtaIndexInSModule(iSM, iMod, iIphi, iIeta);

    LOGF(debug, "Trying to add cell %d \t ieta = %d\t iphi = %d\t amplitude = %1.3f", absId, ieta, iphi, amp);

    // Loop on the nearest cells around, check the highest energy one,
    // and assign its MC label and the time
    float ampMax = 0.f;
    for (int ietai = ieta - 1; ietai <= ieta + 1; ++ietai) {
      for (int iphii = iphi - 1; iphii <= iphi + 1; ++iphii) {

        // Avoid same cell
        if (iphii == iphi && ietai == ieta) {
          continue;
        }

        // Avoid cells out of SM
        if (ietai < 0 || ietai >= NColumns[iSM] || iphii < 0 || iphii >= NRows[iSM]) {
          continue;
        }

        int absIDi = mGeometry->GetAbsCellIdFromCellIndexes(iSM, iphii, ietai);
        // Try to find the cell that will get energy induced
        float ampi = 0.f;
        size_t indexInCells = 0;
        auto itCell = std::find_if((*mCells).begin(), (*mCells).begin() + nCells, [absIDi](const o2::emcal::Cell& cell) {
          return cell.getTower() == absIDi;
        });

        if (itCell != (*mCells).begin() + nCells) {
          // We found a cell, so let's get the amplitude of that cell
          ampi = itCell->getAmplitude();
          if (ampi <= ampMax) {
            continue; // early continue if the new amplitude is not the biggest one
          }
          indexInCells = std::distance((*mCells).begin(), itCell);
          LOGF(debug, "Found cell with index %d", indexInCells);
        } else {
          continue;
        }

        // Remove cells with no energy
        if (ampi <= MinCellEnergy) {
          continue;
        }

        // Only same TCard
        if (!std::get<0>(mGeometry->areAbsIDsFromSameTCard(absId, absIDi))) {
          continue;
        }

        std::vector<int32_t> mclabeli = {(*mCellLabels)[indexInCells].GetLeadingMCLabel()};
        float timei = (*mCells)[indexInCells].getTimeStamp();

        ampMax = ampi;
        mclabel = mclabeli;
        time = timei;
      } // loop phi
    } // loop eta
    // End Assign MC label
    LOGF(debug, "Final ampMax %1.2f\n", ampMax);
    LOGF(debug, "--- End  : Added cell ID %d, ieta %d, iphi %d, E %1.3f, time %1.3e, mc label %d\n", absId, ieta, iphi, amp, time, mclabel[0]);

    // Add the new cell
    (*mCells).emplace_back(absId, amp, time, o2::emcal::intToChannelType(1));
    (*mCellLabels).emplace_back(std::vector<int32_t>{mclabel[0]}, std::vector<float>{0.f});
  } // loop over cells
}

bool EMCCrossTalk::run()
{
  // START PROCESSING
  // Test if cells present
  if ((*mCells).size() == 0) {
    LOGF(error, "No EMCAL cells found, exiting EMCCrossTalk::run()!");
    return false;
  }

  // CELL CROSSTALK EMULATION
  // Compute the induced cell energies by T-Card correlation emulation, ONLY MC
  makeCellTCardCorrelation();

  // Add to existing cells the found induced energies in MakeCellTCardCorrelation() if new signal is larger than 10 MeV.
  addInducedEnergiesToExistingCells();

  // Add new cells with found induced energies in MakeCellTCardCorrelation() if new signal is larger than 10 MeV.
  addInducedEnergiesToNewCells();

  resetArrays();

  return true;
}
