// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \brief  Defines the BaseSelection class for managing and evaluating multiple selections over multiple observables.
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_BASESELECTION_H_
#define PWGCF_FEMTO_CORE_BASESELECTION_H_

#include "PWGCF/Femto/Core/selectionContainer.h"

#include "fairlogger/Logger.h"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{

/// \class BaseSelection
/// \brief Template class for managing selection criteria across multiple observables.
///
/// This class manages an array of SelectionContainer objects, each corresponding to a specific observable.
/// It evaluates which selections are fulfilled, assembles a final bitmask, and tracks required vs. optional cuts.
///
/// \tparam T Type of observable values (mostly floats).
/// \tparam BitmaskType Type used for internal bitmask operations (e.g., uint32_t, uint64_t).
/// \tparam NumObservables Total number of observables handled.
template <typename T, typename BitmaskType, size_t NumObservables>
class BaseSelection
{
 public:
  /// \brief Default constructor.
  BaseSelection() {}

  /// \brief Destructor
  virtual ~BaseSelection() = default;

  /// \brief Add a static-value based selection for a specific observable.
  /// \param selectionValues Vector of threshold values.
  /// \param observableIndex Index of the observable.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold in the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory or optional.
  void addSelection(std::vector<T> const& selectionValues, int observableIndex, limits::LimitType limitType, bool skipMostPermissiveBit, bool isMinimalCut)
  {
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    if (skipMostPermissiveBit) {
      mNSelections += selectionValues.size() - 1;
    } else {
      mNSelections += selectionValues.size();
    }
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    // check if any cut is optional
    if (!isMinimalCut) {
      mHasOptionalSelection = true;
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionValues, limitType, skipMostPermissiveBit, isMinimalCut);
  }

  /// \brief Add a function-based selection for a specific observable.
  /// \param baseName Base name for TF1 functions.
  /// \param lowerLimit Lower bound for the TF1 domain.
  /// \param upperLimit Upper bound for the TF1 domain.
  /// \param selectionValues Function definitions as strings.
  /// \param observableIndex Index of the observable.
  /// \param limitType Type of limit.
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold in the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory or optional.
  void addSelection(std::string const& baseName,
                    T lowerLimit,
                    T upperLimit,
                    std::vector<std::string> const& selectionValues,
                    int observableIndex,
                    limits::LimitType limitType,
                    bool skipMostPermissiveBit,
                    bool isMinimalCut)
  {
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    if (skipMostPermissiveBit) {
      mNSelections += selectionValues.size() - 1;
    } else {
      mNSelections += selectionValues.size();
    }
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(baseName, lowerLimit, upperLimit, selectionValues, limitType, skipMostPermissiveBit, isMinimalCut);
  }

  /// \brief Update the limits of a function-based selection for a specific observable.
  /// \param observable Index of the observable.
  /// \param value Value at which to evaluate the selection functions.
  void updateLimits(int observable, T value) { mSelectionContainers.at(observable).updateLimits(value); }

  /// \brief Reset the internal bitmask and evaluation flags before evaluating a new event.
  void reset()
  {
    mFinalBitmask.reset();
    mPassesMinimalSelections = true;
    // will be true if no optional cut as been defined and
    // will be set to false if we have optional cuts (but will be set to true in the case at least one optional cut succeeds)
    mPassesOptionalSelections = !mHasOptionalSelection;
  }

  /// \brief Evaluate a single observable against its configured selections.
  /// \param observableIndex Index of the observable.
  /// \param value Value of the observable.
  void evaluateObservable(int observableIndex, T value)
  {
    // if there are no selections configured, bail out
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }
    // if any previous observable did not pass minimal selections, there is no point in setting bitmask for other observables
    // minimal selection for each observable is computed after adding it
    if (mPassesMinimalSelections == false) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(value);
    // check if minimal selction for this observable holds
    if (mSelectionContainers.at(observableIndex).passesAsMinimalCut() == false) {
      mPassesMinimalSelections = false;
    }
    // check if any optional selection holds
    if (mSelectionContainers.at(observableIndex).passesAsOptionalCut() == true) {
      mPassesOptionalSelections = true;
    }
  }

  /// \brief Check if all required (minimal) and optional cuts are passed.
  /// \return True if all required and at least one optional cut (if present) is passed.
  bool passesAllRequiredSelections() const
  {
    return mPassesMinimalSelections && mPassesOptionalSelections;
  }

  /// \brief Check if the optional selection for a specific observable is passed.
  /// \param observableIndex Index of the observable.
  /// \return True if at least one optional selection is fulfilled.
  bool passesOptionalSelection(int observableIndex) const
  {
    return mSelectionContainers.at(observableIndex).passesAsOptionalCut();
  }

  /// \brief Assemble the global selection bitmask from individual observable selections.
  void assembleBitmask()
  {
    // if minimal selections are not passed, just set bitmask to 0
    if (mPassesMinimalSelections == false) {
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

  /// \brief Retrieve the assembled bitmask as an integer value.
  /// \return The combined selection bitmask.
  BitmaskType getBitmask() const { return static_cast<BitmaskType>(mFinalBitmask.to_ullong()); }

  /// \brief Retrieve the assembled bitmask as an integer value.
  /// \return The combined selection bitmask.
  BitmaskType getBitmask(int observableIndex) const { return static_cast<BitmaskType>(mSelectionContainers.at(observableIndex).getBitmask().to_ullong()); }

  /// \brief Retrieve the assembled bitmask as an integer value.
  /// \return The combined selection bitmask.
  template <typename R>
  void setBitmask(int observableIndex, R bitmask)
  {
    mSelectionContainers.at(observableIndex).setBitmask(bitmask);
  }

  /// \brief Print detailed information about all configured selections.
  /// \tparam MapType Type used in the observable name map (usually an enum or int).
  /// \param objectName Name of the current object (e.g. particle species).
  /// \param observableNames Map from observable index to human-readable names.
  template <typename MapType>
  void printSelections(const std::string& objectName, const std::unordered_map<MapType, std::string>& observableNames) const
  {
    LOG(info) << "Printing Configuration of " << objectName;

    size_t globalBitIndex = 0; // Tracks bit position across all containers

    for (size_t idx = mSelectionContainers.size(); idx-- > 0;) {
      const auto& container = mSelectionContainers[idx];
      if (container.isEmpty()) {
        continue;
      }

      const MapType key = static_cast<MapType>(idx);
      const std::string& name = observableNames.count(key) ? observableNames.at(key) : "[Unknown]";

      LOG(info) << "Observable: " << name << " (index " << idx << ")";
      LOG(info) << "  Limit type           : " << container.getLimitTypeAsString();
      LOG(info) << "  Minimal cut          : " << (container.isMinimalCut() ? "yes" : "no");
      LOG(info) << "  Skip most permissive : " << (container.skipMostPermissiveBit() ? "yes" : "no");
      LOG(info) << "  Bitmask shift        : " << container.getShift();
      LOG(info) << "  Selections           : ";

      const auto& values = container.getSelectionValues();
      const auto& functions = container.getSelectionFunction();
      const bool useFunctions = !functions.empty();
      const size_t numSelections = useFunctions ? functions.size() : values.size();
      const bool skipMostPermissive = container.skipMostPermissiveBit();

      int valWidth = 20;
      int bitWidth = 30;

      for (size_t j = 0; j < numSelections; ++j) {
        std::stringstream line;

        // Selection string (either value or function)
        const std::string& sel = useFunctions ? std::string(functions[j].GetFormula()->GetExpFormula().Data()) : std::to_string(values[j]);
        line << "    " << std::left << std::setw(valWidth) << sel;

        // Bitmask
        if (skipMostPermissive && j == 0) {
          line << std::setw(bitWidth) << "-> loosest minimal selection, no bit saved";
        } else {
          const uint64_t bitmask = uint64_t{1} << globalBitIndex++;
          line << std::setw(bitWidth) << ("-> bitmask: " + std::to_string(bitmask));
        }

        LOG(info) << line.str();
      }

      LOG(info) << ""; // blank line between observables
    }
    LOG(info) << "Printing done";
  }

 protected:
  std::array<SelectionContainer<T, BitmaskType>, NumObservables> mSelectionContainers = {}; ///< Array containing all selections
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mFinalBitmask = {};                           ///< final bitmaks
  size_t mNSelections = 0;                                                                  ///< Number of selections
  bool mPassesMinimalSelections = true;                                                     ///< Set to true if all minimal (mandatory) selections are passed
  bool mHasOptionalSelection = false;                                                       ///< Set to true if at least one selections is optional
  bool mPassesOptionalSelections = true;                                                    ///< Set to true if at least one optional (non-mandatory) selections is passed
};
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_BASESELECTION_H_
