// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file baseSelection.h
/// \brief Definition of the BaseSelection class
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_

#include "PWGCF/FemtoUnited/Core/selectionContainer.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>

namespace o2::analysis::femtounited
{
/// \class BaseSelection
/// \brief Base class to contain all cuts and assemble bitmask
/// \tparam T Data type used for the selections (float/int/...)
/// \tparam BitmaskType Compute size of the bitmask from BitmaskType
/// \tparam NObservables Number of observables
template <typename T, typename BitmaskType, size_t NumObservables>
class BaseSelection
{
 public:
  /// Constructor
  BaseSelection() {}

  /// Destructor
  virtual ~BaseSelection() = default;

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \param configSelections Vector of configurables containing the values employed for the selection
  /// \param Observable Observable to be employed for the selection
  /// \param limitType Type of the selection limit
  void addSelection(std::vector<T> const& selectionValues, int observableIndex, limits::LimitType limitType, bool skipMostPermissiveBit, bool isMinimalCut)
  {
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    mNSelections += selectionValues.size();
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    // check if any cut is optional
    if (!isMinimalCut) {
      mHasOptionalCuts = true;
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionValues, limitType, skipMostPermissiveBit, isMinimalCut);
  }

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \param configSelection Vector from configurable containing the values employed for the selection
  /// \param observableType Observable to be employed for the selection
  /// \param limitType Type of the selection limit
  void addSelection(std::string const& baseName, T lowerLimit, T upperLimit, std::vector<std::string> const& selectionValues, int observableIndex, limits::LimitType limitType, bool skipMostPermissiveBit, bool isMinimalCut)
  {
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    mNSelections += selectionValues.size();
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(baseName, lowerLimit, upperLimit, selectionValues, limitType, skipMostPermissiveBit, isMinimalCut);
  }

  void updateLimits(int observable, T value) { mSelectionContainers.at(observable).updateLimits(value); }

  void reset()
  {
    mFinalBitmask.reset();
    mPassesMinimalCuts = true;
    // will be true if no optional cut as been defined and
    // will be set to false if we have optional cuts (but will be set to true in the case at least one optional cut succeeds)
    mPassesOptionalCuts = !mHasOptionalCuts;
  }

  /// set bitmask for a given observable
  /// \param observable Observable to be checked
  /// \param value Value of the observable
  void evaluateObservable(int observableIndex, T value)
  {
    // if there are no selections configured, bail out
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }
    // if any previous observable did not pass minimal selections, there is no point in setting bitmask for other observables
    // minimal selection for each observable is computed after adding it
    if (mPassesMinimalCuts == false) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(value);
    // check if minimal selction for this observable holds
    if (mSelectionContainers.at(observableIndex).passesAsMinimalCut() == false) {
      mPassesMinimalCuts = false;
    }
    // check if any optional selection holds
    if (mSelectionContainers.at(observableIndex).passesAsOptionalCut() == true) {
      mPassesOptionalCuts = true;
    }
  }

  /// check if minimal Selections are passed
  bool passesAllRequiredSelections() const
  {
    return mPassesMinimalCuts && mPassesOptionalCuts;
  }

  bool passesOptionalCut(int observableIndex) const
  {
    return mSelectionContainers.at(observableIndex).passesAsOptionalCut();
  }

  /// assemble final bitmask
  void assembleBitmask()
  {
    // if minimal selections are not passed, just set bitmask to 0
    if (mPassesMinimalCuts == false) {
      mFinalBitmask.reset();
      return;
    }

    // to assemble bitmask, convert all bitmask into integers
    // shift the current one and add the new bits
    for (auto const& selectionContainer : mSelectionContainers) {
      // if there are no selections for a certain observable, skip
      if (selectionContainer.isEmpty()) {
        continue;
      }
      // Shift the result to make space and add the new value
      mFinalBitmask = (mFinalBitmask << selectionContainer.getShift()) | selectionContainer.getBitmask();
    }
  }

  BitmaskType getBitmask() const { return static_cast<BitmaskType>(mFinalBitmask.to_ullong()); }

#include <iomanip> // for setw, left
#include <sstream>

  template <typename MapType>
  void printSelections(const std::string& objectName, const std::unordered_map<MapType, std::string>& observableNames) const
  {
    LOG(info) << "**************************************** FemtoProducer ****************************************";
    LOG(info) << objectName << "\n";

    size_t globalBitIndex = 0; // Track absolute bit position across all containers

    for (size_t idx = mSelectionContainers.size(); idx-- > 0;) {
      const auto& container = mSelectionContainers[idx];
      if (container.isEmpty()) {
        continue;
      }

      std::string name = "[Unknown]";
      auto key = static_cast<MapType>(idx);
      if (observableNames.count(key)) {
        name = observableNames.at(key);
      }

      LOG(info) << "Observable: " << name << " (index " << idx << ")";
      LOG(info) << "  Limit type           : " << container.getLimitTypeAsString();
      LOG(info) << "  Values (with bit indices and bitmask):";

      const auto& values = container.getSelectionValues();
      bool skipMostPermissive = container.skipMostPermissiveBit();

      constexpr int valWidth = 15;
      constexpr int bitWidth = 20;
      constexpr int maskWidth = 12;

      for (size_t j = 0; j < values.size(); ++j) {
        std::stringstream line;
        line << "    "
             << std::left << std::setw(valWidth) << values[j];

        if (skipMostPermissive && j == 0) {
          line << std::setw(bitWidth) << "[no bit (skipped)]"
               << std::setw(maskWidth) << "";
        } else {
          uint64_t bitmask = uint64_t{1} << globalBitIndex;
          line << std::setw(bitWidth) << ("[bit " + std::to_string(globalBitIndex) + "]")
               << " bitmask: " << bitmask;
          ++globalBitIndex;
        }

        LOG(info) << line.str();
      }

      LOG(info) << "  Minimal cut          : " << (container.isMinimalCut() ? "yes" : "no");
      LOG(info) << "  Skip most permissive : " << (skipMostPermissive ? "yes" : "no");
      LOG(info) << "  Bitmask shift        : " << container.getShift();
      LOG(info) << ""; // blank line between observables
    }

    LOG(info) << "***********************************************************************************************";
  }

 protected:
  std::array<SelectionContainer<T, BitmaskType>, NumObservables> mSelectionContainers = {}; ///< Array containing all selections
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mFinalBitmask = {};                           ///< final bitmaks
  size_t mNSelections = 0;                                                                  ///< Number of selections
  bool mPassesMinimalCuts = true;                                                           ///< Set to true if all minimal (mandatory) cuts are passed
  bool mHasOptionalCuts = false;                                                            ///< Set to true if at least one cut is optional
  bool mPassesOptionalCuts = true;                                                          ///< Set to true if at least one optional (non-mandatory) cut is passed
};
} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
