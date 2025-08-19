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

/// \file emcalCrossTalkEmulation.h
/// \brief emulation of emcal cross talk for simulations
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe-University

#ifndef PWGJE_CORE_EMCALCROSSTALKEMULATION_H_
#define PWGJE_CORE_EMCALCROSSTALKEMULATION_H_

#include <DataFormatsEMCAL/Cell.h>
#include <DataFormatsEMCAL/CellLabel.h>
#include <EMCALBase/Geometry.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>

#include <TH1.h>
#include <TRandom3.h>

#include <array>
#include <string>
#include <vector>

namespace o2::emccrosstalk
{
// cell types for enegery induction
enum InductionCellType {
  UpDown = 0,
  UpDownLeftRight,
  LeftRight,
  Up2Down2,
  NInductionCellType
};

// default values for cross talk emulation
// small value to use for comarison equal to 0, std::abs(x) > epsilon
static constexpr float Epsilon = 1e-6f;

// default for inducedEnergyLossConstant
// static constexpr float DefaultIELC[1][4] = {{0.02f, 0.02f, 0.02f, 0.f}}; // default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml (but is deactivated per default)
static constexpr float DefaultIELC[1][4] = {{0.f, 0.f, 0.f, 0.f}}; // default from pp 13 TeV

// default for inducedEnergyLossFraction, default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml same as in https://cds.cern.ch/record/2910556/
static constexpr float DefaultIELF[20][4] = {{1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.20e-02, 1.20e-02, 1.20e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.20e-02, 1.20e-02, 1.20e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.20e-02, 1.20e-02, 1.20e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {1.20e-02, 1.20e-02, 1.20e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {1.15e-02, 1.15e-02, 1.15e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f},
                                             {0.80e-02, 0.80e-02, 0.80e-02, 0.f}};

// default for inducedEnergyLossFractionP1, default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml same as in https://cds.cern.ch/record/2910556/
static constexpr float DefaultIELFP1[1][4] = {{-1.1e-03, -1.1e-03, -1.1e-03, 0.f}};

// default for inducedEnergyLossFractionWidth, default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml same as in https://cds.cern.ch/record/2910556/
static constexpr float DefaultIELFWidth[1][4] = {{5.0e-03, 5.0e-03, 5.0e-03, 0.f}};

// default for inducedEnergyLossMinimumFractionCentralEta, IF someone wants to try it. This is purely to document those test numbers from AliEmcalCorrectionConfiguration.yaml! Default is to not use this! Values from default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml
// std::vector<float> DefaultIELMFCE = {6.8e-3, 7.5e-3, 6.8e-3, 9.0e-3, 6.8e-3, 6.8e-3, 6.8e-3, 9.0e-3, 5.2e-3, 5.2e-3, 7.5e-3, 6.8e-3, 6.8e-3, 6.8e-3, 5.2e-3, 5.2e-3, 6.8e-3, 5.2e-3, 5.2e-3, 5.2e-3};

struct EmcCrossTalkConf : o2::framework::ConfigurableGroup {
  std::string prefix = "emccrosstalk";
  o2::framework::Configurable<bool> enableCrossTalk{"enableCrossTalk", false, "Flag to enable cross talk emulation. This should only ever be used for MC!"};
  o2::framework::Configurable<bool> createHistograms{"createHistograms", false, "Flag to enable QA histograms."};
  o2::framework::Configurable<bool> printConfiguration{"printConfiguration", true, "Flag to print the configuration after initialization."};
  o2::framework::Configurable<bool> conserveEnergy{"conserveEnergy", true, "Flag to enable cluster energy conservation."};
  o2::framework::Configurable<bool> randomizeTCardInducedEnergy{"randomizeTCardInducedEnergy", true, "Flag to randomize the energy fraction induced by the TCard."};
  o2::framework::Configurable<float> inducedTCardMinimumCellEnergy{"inducedTCardMinimumCellEnergy", 0.01f, "Minimum cell energy in GeV induced by the TCard."};
  o2::framework::Configurable<float> inducedTCardMaximum{"inducedTCardMaximum", 100.f, "Maximum energy in GeV induced by the TCard."};
  o2::framework::Configurable<float> inducedTCardMinimum{"inducedTCardMinimum", 0.1f, "Minimum energy in GeV induced by the TCard + cell energy, IMPORTANT use the same value as the clusterization cell E threshold or not too far from it."};
  o2::framework::Configurable<float> inducedTCardMaximumELeak{"inducedTCardMaximumELeak", 0.f, "Maximum energy in GeV that is going to be leaked independently of what is set with inducedTCardMinimum."};
  // For each of the following Array2D configurables we will use:
  // empty vector == disabled
  // vector of size 4 == same for all SM
  // vector of vectors with size nSM * 4 == each SM has its own setting
  // the 4 values (0-3) correspond to in relative [row,col]: 0: [+-1,0], 1: [+-1,+or-1], 2:  [0,+or-1], 3: [+-2, 0 AND +or-1]
  // +---+---+-----+---+---+--+--+--+
  // | 3 | 0 | Hit | 0 | 3 |  |  |  |
  // +---+---+-----+---+---+--+--+--+
  // | 3 | 1 |  2  | 1 | 3 |  |  |  |
  // +---+---+-----+---+---+--+--+--+
  // For the std::vector<float> it is similar, empty vector means not used, single value means one value for all SM and 20 values means specifiyng a value for all SM
  o2::framework::Configurable<o2::framework::Array2D<float>> inducedEnergyLossConstant{"inducedEnergyLossConstant", {DefaultIELC[0], 1, 4}, "Constant energy lost by max energy cell in one of T-Card cells. Empty vector == disabled, size 4 vector == enabled. For information on the exact formatting please check the header file."};
  o2::framework::Configurable<o2::framework::Array2D<float>> inducedEnergyLossFraction{"inducedEnergyLossFraction", {DefaultIELF[0], 20, 4}, "Fraction of energy lost by max energy cell in one of T-Card cells."};
  o2::framework::Configurable<o2::framework::Array2D<float>> inducedEnergyLossFractionP1{"inducedEnergyLossFractionP1", {DefaultIELFP1[0], 1, 4}, "Slope parameter of fraction of energy lost by max energy cell in one of T-Card cells."};
  o2::framework::Configurable<o2::framework::Array2D<float>> inducedEnergyLossFractionWidth{"inducedEnergyLossFractionWidth", {DefaultIELFWidth[0], 1, 4}, "Fraction of energy lost by max energy cell in one of T-Card cells, width of random gaussian."};
  // default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml :
  // o2::framework::Configurable<std::vector<float>> inducedEnergyLossMinimumFraction{"inducedEnergyLossMinimumFraction", {3.5e-3f, 5.0e-3f, 4.5e-3f, 6.0e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 6.0e-3f, 3.5e-3f, 3.5e-3f, 5.0e-3f, 5.0e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f, 3.5e-3f}, "Minimum induced energy fraction when linear dependency is set."};
  // value from https://cds.cern.ch/record/2910556/:
  o2::framework::Configurable<std::vector<float>> inducedEnergyLossMinimumFraction{"inducedEnergyLossMinimumFraction", {2.35e-3f, 2.5e-3f, 2.35e-3f, 3.0e-3f, 2.35e-3f, 2.35e-3f, 2.35e-3f, 3.0e-3f, 1.75e-3f, 1.75e-3f, 2.5e-3f, 2.35e-3f, 2.35e-3f, 2.35e-3f, 1.75e-3f, 1.75e-3f, 2.35e-3f, 1.75e-3f, 1.75e-3f, 1.75e-3f}, "Minimum induced energy fraction when linear dependency is set."};

  o2::framework::Configurable<std::vector<float>> inducedEnergyLossMinimumFractionCentralEta{"inducedEnergyLossMinimumFractionCentralEta", {}, "Minimum induced energy fraction when linear dependency is set. For |eta| < 0.22, if empty no difference in eta. NOT TUNED for TESTING!"};

  // default from https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml :
  // o2::framework::Configurable<std::vector<float>> inducedEnergyLossMaximumFraction{"inducedEnergyLossMaximumFraction", {0.018f}, "Maximum induced energy fraction when linear dependency is set."};
  // value from https://cds.cern.ch/record/2910556/:
  o2::framework::Configurable<std::vector<float>> inducedEnergyLossMaximumFraction{"inducedEnergyLossMaximumFraction", {0.016f, 0.016f, 0.016f, 0.018f, 0.016f, 0.016f, 0.016f, 0.018f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f, 0.016f}, "Maximum induced energy fraction when linear dependency is set."};
  o2::framework::Configurable<std::vector<float>> inducedEnergyLossProbability{"inducedEnergyLossProbability", {1.0f}, "Fraction of times max cell energy correlates with cross cells."};
};

static constexpr int NSM = 20;                // Number of Supermodules (12 for EMCal + 8 for DCal)
static constexpr int NCells = 17664;          // Number of cells in the EMCal
static constexpr int NNeighbourCases = 4;     // 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells
static constexpr int FirstDCal23SM = 12;      // index of the first 2/3 DCal SM
static constexpr int LastDCal23SM = 17;       // index of the last 2/3 DCal SM
static constexpr float MinCellEnergy = 0.01f; // Minimum energy a new cell needs to be added

static constexpr int NColumns[NSM] = {48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 32, 32, 32, 32, 32, 32, 48, 48};
static constexpr int NRows[NSM] = {24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 8, 8, 24, 24, 24, 24, 24, 24, 8, 8};

// these labels are for later once labeledArrays work on hyperloop. Currently they sadly only allow fixed size not variable size.
// static const std::vector<std::string> labelsSM{"SM0/all", "SM1", "SM2", "SM3", "SM4", "SM5", "SM6", "SM7", "SM8", "SM9", "SM10", "SM11", "SM12", "SM13", "SM14", "SM15", "SM16", "SM17", "SM18", "SM19"};
// static const std::vector<std::string> labelsCells = {"Up&Down", "Up&Down x Left|Right", "Left|Right", "2Up&Down + 2Up&Down xLeft|Right"};

class EMCCrossTalk
{

 public:
  ~EMCCrossTalk()
  {
    LOG(info) << "Destroying EMCCrossTalk";
  }

  /// \brief Basic init function.
  /// \param config configurable group containing the config for the cross talk emulation
  void initObjects(const EmcCrossTalkConf& config);

  /// \brief Reset arrays containing information for all possible cells.
  /// \details mTCardCorrCellsEner and mTCardCorrCellsNew
  void resetArrays();

  /// \brief Sets the pointer the current vector of cells.
  /// \param cells pointer to emcal cells of the current event
  /// \param cellLabels pointer to emcal cell labels of the current event
  void setCells(std::vector<o2::emcal::Cell>& cells, std::vector<o2::emcal::CellLabel>& cellLabels);

  /// \brief Main function to call later to perform the full cross talk emulation
  /// \return flag if everything went well or not
  bool run();

  /// \brief Recover each cell amplitude and absId and induce energy in cells in cross of the same T-Card
  void makeCellTCardCorrelation();

  /// \brief Add to existing cells the found induced energies in makeCellTCardCorrelation() if new signal is larger than 10 MeV.
  /// \details Need to destroy/create the default cells list and do a copy from the old to the new via a temporal array fAODCellsTmp. Not too nice or fast, but it works.
  void addInducedEnergiesToExistingCells();

  /// \brief Add new cells with found induced energies in makeCellTCardCorrelation() if new signal is larger than 10 MeV.
  void addInducedEnergiesToNewCells();

  /// \brief Calculate the induced energy in a cell belonging to thesame T-Card as the reference cell. Used in makeCellTCardCorrelation()
  /// \param absId Id number of cell in same T-Card as reference cell
  /// \param absIdRef Id number of reference cell
  /// \param iSM Supermodule number of cell
  /// \param ampRef Amplitude of the reference cell
  /// \param cellCase Type of cell with respect reference cell 0: up or down, 1: up or down on the diagonal, 2: left or right, 3: 2nd row up/down both left/right
  void calculateInducedEnergyInTCardCell(int absId, int absIdRef, int iSM, float ampRef, int cellCase);

 private:
  // T-Card correlation emulation, do on MC
  bool mTCardCorrClusEnerConserv;                // When making correlation, subtract from the reference cell the induced energy on the neighbour cells
  std::array<float, NCells> mTCardCorrCellsEner; //  Array with induced cell energy in T-Card neighbour cells
  std::array<bool, NCells> mTCardCorrCellsNew;   //  Array with induced cell energy in T-Card neighbour cells, that before had no signal

  o2::framework::Array2D<float> mTCardCorrInduceEner;           // Induced energy loss gauss constant on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  o2::framework::Array2D<float> mTCardCorrInduceEnerFrac;       // Induced energy loss gauss fraction param0 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  o2::framework::Array2D<float> mTCardCorrInduceEnerFracP1;     // Induced energy loss gauss fraction param1 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param1
  o2::framework::Array2D<float> mTCardCorrInduceEnerFracWidth;  // Induced energy loss gauss witdth on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells
  std::array<float, NSM> mTCardCorrInduceEnerFracMax;           // In case fTCardCorrInduceEnerFracP1  is non null, restrict the maximum fraction of induced energy per SM
  std::array<float, NSM> mTCardCorrInduceEnerFracMin;           // In case fTCardCorrInduceEnerFracP1  is non null, restrict the minimum fraction of induced energy per SM
  std::array<float, NSM> mTCardCorrInduceEnerFracMinCentralEta; // In case fTCardCorrInduceEnerFracP1  is non null, restrict the minimum fraction of induced energy per SM. Different at central |eta| < 0.22
  std::array<float, NSM> mTCardCorrInduceEnerProb;              // Probability to induce energy loss per SM

  TRandom3 mRandom;     //  Random generator
  bool mRandomizeTCard; //  Use random induced energy

  float mTCardCorrMinAmp;          //  Minimum cell energy to induce signal on adjacent cells
  float mTCardCorrMinInduced;      //  Minimum induced energy signal on adjacent cells, sum of induced plus original energy, use same as cell energy clusterization cut
  float mTCardCorrMaxInducedELeak; //  Maximum value of induced energy signal that is always leaked, ~5-10 MeV
  float mTCardCorrMaxInduced;      //  Maximum induced energy signal on adjacent cells

  std::vector<o2::emcal::Cell>* mCells = nullptr;           // Pointer to the original cells of the current event
  std::vector<o2::emcal::CellLabel>* mCellLabels = nullptr; // Pointer to the original cell labels of the current event
  std::vector<o2::emcal::Cell> mCellsTmp;                   // Temporal vector of cells (copy)

  o2::emcal::Geometry* mGeometry; // EMCal geometry
};

} // namespace o2::emccrosstalk

#endif // PWGJE_CORE_EMCALCROSSTALKEMULATION_H_
