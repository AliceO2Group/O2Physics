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

#include "Framework/HistogramRegistry.h"

#include "fairlogger/Logger.h"

#include <climits>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>
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
  BaseSelection() = default;

  /// \brief Destructor
  virtual ~BaseSelection() = default;

  /// \brief Add a static-value based selection for a specific observable.
  /// \param selectionValues Vector of threshold values.
  /// \param observableIndex Index of the observable.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold in the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory or optional.
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    std::vector<T> const& selectionValues,
                    limits::LimitType limitType,
                    bool skipMostPermissiveBit,
                    bool isMinimalCut,
                    bool isOptionCut)
  {
    // check index
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    // init selection container for selection at given index
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, selectionValues, limitType, skipMostPermissiveBit, isMinimalCut, isOptionCut);

    // check if any selections are configured
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }

    // keep track of selections and bits
    mNSelectionBits += mSelectionContainers.at(observableIndex).getShift();
    mNSelection += mSelectionContainers.at(observableIndex).getNSelections();

    if (mNSelectionBits > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " number of bits are supported";
    }
    // check if any selection is minimal
    if (mSelectionContainers.at(observableIndex).isMinimalCut()) {
      mHasMinimalSelection = true;
    }
    // check if selection is optional
    if (mSelectionContainers.at(observableIndex).isOptionalCut()) {
      mHasOptionalSelection = true;
    }
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
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    T lowerLimit,
                    T upperLimit,
                    std::vector<std::string> const& functions,
                    limits::LimitType limitType,
                    bool skipMostPermissiveBit,
                    bool isMinimalCut,
                    bool isOptionalCut)
  {
    if (static_cast<size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, lowerLimit, upperLimit, functions, limitType, skipMostPermissiveBit, isMinimalCut, isOptionalCut);

    // check if any selections are configured
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }

    // advance mNSelections so we can use it as offset for next selection
    mNSelectionBits += mSelectionContainers.at(observableIndex).getShift();
    mNSelection += mSelectionContainers.at(observableIndex).getNSelections();

    if (mNSelectionBits > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    // keep track of selection selections
    // check if any cut is minimal
    if (mSelectionContainers.at(observableIndex).isMinimalCut()) {
      mHasMinimalSelection = true;
    }
    // check if any selection is optional
    if (mSelectionContainers.at(observableIndex).isOptionalCut()) {
      mHasOptionalSelection = true;
    }
  }

  /// \brief Add a boolean based selection for a specific observable.
  /// \param mode Whether the selection is not applied, minimal or optional cut
  /// \param observableIndex Index of the observable.
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    int mode)
  {
    switch (mode) {
      case -1: // cut is optional and we store bit for the cut
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{1}, limits::LimitType::kEqual, false, false, true);
        mHasOptionalSelection = true;
        mNSelectionBits += 1;
        mNSelection += 1;
        break;
      case 0: // cut is not applied, initalize with empty vector, so we bail out later
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{}, limits::LimitType::kEqual, false, false, false);
        break;
      case 1: // cut is added as mininal selection (since it is only one value, no extra bit is stored)
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{1}, limits::LimitType::kEqual, true, true, false);
        mHasMinimalSelection = true;
        mNSelection += 1;
        break;
      default:
        LOG(fatal) << "Invalid switch for boolean selection";
    }
    if (mNSelectionBits > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
  }

  /// \brief Update the limits of a function-based selection for a specific observable.
  /// \param observable Index of the observable.
  /// \param value Value at which to evaluate the selection functions.
  void updateLimits(int observable, T value)
  {
    mSelectionContainers.at(observable).updateLimits(value);
  }

  /// \brief Reset the internal bitmask and evaluation flags before evaluating a new event.
  void reset()
  {
    mFinalBitmask.reset();
    for (std::size_t i = 0; i < mSelectionContainers.size(); i++) {
      mSelectionContainers.at(i).reset();
    }
    if (mHasMinimalSelection) {
      mPassesMinimalSelections = true;
    }
    if (mHasOptionalSelection) {
      mPassesOptionalSelections = false;
    }
  }

  void reset(int observableIndex) { mSelectionContainers.at(observableIndex).reset(); }

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
    if (!mPassesMinimalSelections) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(value);
    // check if minimal selction for this observable holds
    // if one minimal selection is not fullfilled, the condition failes
    if (mHasMinimalSelection) {
      if (!mSelectionContainers.at(observableIndex).passesAsMinimalCut()) {
        mPassesMinimalSelections = false;
      }
    }
    // check if any optional selection holds
    // if one optional selection is fullfilled, the condition succeeds
    if (mHasOptionalSelection) {
      if (mSelectionContainers.at(observableIndex).passesAsOptionalCut()) {
        mPassesOptionalSelections = true;
      }
    }
  }

  /// \brief Evaluate a single observable against its configured selections.
  /// \param observableIndex Index of the observable.
  /// \param values vector of values of the observable.
  void evaluateObservable(int observableIndex, std::vector<T> values)
  {
    // if there are no selections configured, bail out
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }
    // if any previous observable did not pass minimal selections, there is no point in setting bitmask for other observables
    // minimal selection for each observable is computed after adding it
    if (!mPassesMinimalSelections) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(values);
    // check if minimal selction for this observable holds
    if (mHasMinimalSelection) {
      if (mSelectionContainers.at(observableIndex).passesAsMinimalCut() == false) {
        mPassesMinimalSelections = false;
      }
    }
    // check if any optional selection holds
    if (mHasOptionalSelection) {
      if (mSelectionContainers.at(observableIndex).passesAsOptionalCut() == true) {
        mPassesOptionalSelections = true;
      }
    }
  }

  /// \brief Add comments to specific observabel
  void addComments(int observableIndex, std::vector<std::string> const& comments) { mSelectionContainers.at(observableIndex).addComments(comments); }

  /// \brief Check if all required (minimal) and optional cuts are passed.
  /// \return True if all required and at least one optional cut (if present) is passed.
  bool passesAllRequiredSelections() const
  {
    if (mHasMinimalSelection && !mHasOptionalSelection) {
      return mPassesMinimalSelections;
    }
    if (!mHasMinimalSelection && mHasOptionalSelection) {
      return mPassesOptionalSelections;
    }
    if (mHasMinimalSelection && mHasOptionalSelection) {
      return mPassesMinimalSelections && mPassesOptionalSelections;
    }
    return true;
  }

  /// \brief Check if the optional selection for a specific observable is passed.
  /// \param observableIndex Index of the observable.
  /// \return True if at least one optional selection is fulfilled.
  bool passesOptionalSelection(int observableIndex) const
  {
    return mSelectionContainers.at(observableIndex).passesAsOptionalCut();
  }

  /// \brief Assemble the global selection bitmask from individual observable selections.
  template <const char* HistName>
  void assembleBitmask()
  {
    mHistRegistry->fill(HIST(HistName), mNSelection);
    // if the required selections are not passed, we can break early
    if (!this->passesAllRequiredSelections()) {
      mFinalBitmask.reset();
      return;
    }
    mHistRegistry->fill(HIST(HistName), mNSelection + 1);

    int binCenter = 0;
    // to assemble bitmask, convert all bitmask into integers
    // shift the current one and add the new bits
    for (auto const& selectionContainer : mSelectionContainers) {
      // if there are no selections for a certain observable, skip
      if (selectionContainer.isEmpty()) {
        continue;
      }
      // Shift the result to its offset and add the new values
      mFinalBitmask |= (selectionContainer.getBitmask() << selectionContainer.getOffset());

      for (int j = 0; j < selectionContainer.getNSelections(); ++j) {
        if (j == 0 && selectionContainer.isMinimalCut()) {
          // minimal cuts are always filled
          mHistRegistry->fill(HIST(HistName), binCenter);
        } else {
          // use container's internal offset for checking the bit
          if (mFinalBitmask.test(selectionContainer.getBitPosition(j))) {
            mHistRegistry->fill(HIST(HistName), binCenter);
          }
        }
        binCenter++;
      }
    }
  }

  /// \brief Retrieve the assembled bitmask as an integer value.
  /// \return The combined selection bitmask.
  BitmaskType getBitmask() const { return static_cast<BitmaskType>(mFinalBitmask.to_ullong()); }

  /// \brief Retrieve the assembled bitmask as an integer value.
  /// \return The combined selection bitmask.
  BitmaskType getBitmask(int observableIndex) const { return static_cast<BitmaskType>(mSelectionContainers.at(observableIndex).getBitmask().to_ullong()); }

  /// \brief Set the assembled bitmask for on observable
  /// \return The combined selection bitmask.
  template <typename R>
  void setBitmask(int observableIndex, R bitmask)
  {
    mSelectionContainers.at(observableIndex).setBitmask(bitmask);
  }

  T getLoosestSelection(int observableIndex) const { return mSelectionContainers.at(observableIndex).getLoosestSelection(); }

  void printSelections(const std::string& objectName) const
  {
    LOG(info) << "Printing Configuration of " << objectName;
    for (size_t idx = 0; idx < mSelectionContainers.size(); ++idx) {
      const auto& container = mSelectionContainers[idx];
      if (container.isEmpty()) {
        continue;
      }

      LOG(info) << "  Observable: " << container.getSelectionName() << " (index " << idx << ")";
      LOG(info) << "  Limit type               : " << container.getLimitTypeAsString();
      LOG(info) << "  Skip most permissive Bit : " << (container.skipMostPermissiveBit() ? "yes" : "no");
      LOG(info) << "  Minimal cut              : " << (container.isMinimalCut() ? "yes" : "no");
      LOG(info) << "  Optional cut             : " << (container.isOptionalCut() ? "yes" : "no");
      LOG(info) << "  Bitmask offset           : " << container.getOffset();
      LOG(info) << "  Bitmask shift            : " << container.getShift();
      LOG(info) << "  Selections:";

      const bool useFunctions = container.isUsingFunctions();
      const auto& values = container.getSelectionValues();
      const auto& functions = container.getSelectionFunction();
      const auto& comments = container.getComments();

      for (int j = 0; j < container.getNSelections(); ++j) {

        std::stringstream line;
        std::string sel = useFunctions ? std::string(functions[j].GetExpFormula().Data()) : std::to_string(values[j]);

        line << "    " << std::left << std::setw(25) << sel;

        if (j == 0 && container.isMinimalCut()) {
          line << "-> minimal cut, no bit saved";
        } else {
          int bit = container.getOffset() + (j - (container.isMinimalCut() ? 1 : 0));
          line << "-> Bit: 0x" << std::hex << std::uppercase << (1ULL << bit) << std::dec;
        }

        if (!comments.empty()) {
          line << " (" << comments.at(j) << ")";
        }
        LOG(info) << line.str();
      }
      LOG(info) << "";
    }
    LOG(info) << "Number of occupied bits: " << mNSelectionBits << " / " << sizeof(BitmaskType) * CHAR_BIT;
    LOG(info) << "Printing done";
  }

  template <const char* HistName>
  void setupContainers(o2::framework::HistogramRegistry* registry)
  {
    mHistRegistry = registry;
    // Create histogram with correct number of bins
    int nBins = mNSelection + 2;
    mHistRegistry->add(HistName, "; Selection Bits; Entries", o2::framework::kTH1F, {{nBins, -0.5, nBins - 0.5}});

    size_t binIndex = 0;
    int offset = 0;
    for (size_t idx = 0; idx < mSelectionContainers.size(); ++idx) {
      auto& container = mSelectionContainers[idx];
      if (container.isEmpty()) {
        continue;
      }
      container.setOffset(offset);
      offset += container.getShift();
      for (int j = 0; j < container.getNSelections(); j++) {
        std::string label = container.getBinLabel(j);
        mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(binIndex + 1, label.c_str());
        binIndex++;
      }
    }
    mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(mNSelection + 1, "All analyzed");
    mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(mNSelection + 2, "All passed");
  }

 protected:
  o2::framework::HistogramRegistry* mHistRegistry = nullptr;
  std::array<SelectionContainer<T, BitmaskType>, NumObservables> mSelectionContainers = {}; ///< Array containing all selections
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mFinalBitmask = {};                           ///< final bitmaks
  std::size_t mNSelectionBits = 0;                                                          ///< Number of selections (all - minimal selections)
  int mNSelection = 0;                                                                      ///< Number of selections all selections
  bool mHasMinimalSelection = false;                                                        ///< Set to true if all minimal (mandatory) selections are passed
  bool mPassesMinimalSelections = true;                                                     ///< Set to true if all minimal (mandatory) selections are passed
  bool mHasOptionalSelection = false;                                                       ///< Set to true if at least one selections is optional
  bool mPassesOptionalSelections = false;                                                   ///< Set to true if at least one optional (non-mandatory) selections is passed
};
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_BASESELECTION_H_
