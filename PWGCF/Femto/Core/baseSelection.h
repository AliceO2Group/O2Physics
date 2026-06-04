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

#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>

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
/// \tparam BitmaskType Integer type used for bitmask operations (e.g., uint32_t, uint64_t).
/// \tparam NumObservables Total number of observables handled.
template <typename T, typename BitmaskType, std::size_t NumObservables>
class BaseSelection
{
 public:
  /// \brief Default constructor.
  BaseSelection() = default;

  /// \brief Destructor
  virtual ~BaseSelection() = default;

  /// \brief Add a static-value based selection for a specific observable.
  /// \param observableIndex Index of the observable.
  /// \param selectionName Name of the selection.
  /// \param selectionValues Vector of threshold values.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory (must be passed for the candidate to be accepted).
  /// \param isOptionalCut Whether this cut is optional (candidate is accepted if any optional cut passes).
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    std::vector<T> const& selectionValues,
                    limits::LimitType limitType,
                    bool skipMostPermissiveBit,
                    bool isMinimalCut,
                    bool isOptionalCut)
  {
    // check index
    if (static_cast<std::size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    // init selection container for selection at given index
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, selectionValues, limitType, skipMostPermissiveBit, isMinimalCut, isOptionalCut);

    init(observableIndex);
  }

  /// \brief Add a function-based selection for a specific observable.
  /// \param observableIndex Index of the observable.
  /// \param selectionName Name of the selection.
  /// \param lowerLimit Lower bound of the TF1 domain.
  /// \param upperLimit Upper bound of the TF1 domain.
  /// \param functions Selection threshold functions as strings (parsed as TF1 expressions).
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory (must be passed for the candidate to be accepted).
  /// \param isOptionalCut Whether this cut is optional (candidate is accepted if any optional cut passes).
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
    if (static_cast<std::size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, lowerLimit, upperLimit, functions, limitType, skipMostPermissiveBit, isMinimalCut, isOptionalCut);

    init(observableIndex);
  }

  /// \brief Add a static-value based selection for a specific observable.
  /// \param observableIndex Index of the observable.
  /// \param selectionName Name of the selection.
  /// \param selectionValues Vector of threshold values.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the loosest threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this cut is mandatory (must be passed for the candidate to be accepted).
  /// \param isOptionalCut Whether this cut is optional (candidate is accepted if any optional cut passes).
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    std::vector<std::string> const& selectionRanges,
                    bool skipMostPermissiveBit,
                    bool isMinimalCut,
                    bool isOptionalCut)
  {
    // check index
    if (static_cast<std::size_t>(observableIndex) >= NumObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NumObservables;
    }
    // init selection container for selection at given index
    mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, selectionRanges, skipMostPermissiveBit, isMinimalCut, isOptionalCut);

    init(observableIndex);
  }

  /// \brief Add a boolean-based selection for a specific observable.
  /// \param observableIndex Index of the observable.
  /// \param selectionName Name of the selection.
  /// \param mode Controls how the cut is applied:
  ///             -1 = optional cut, bit is stored in bitmask;
  ///              0 = cut is disabled, no bit stored;
  ///              1 = minimal (mandatory) cut, no extra bit stored since only one threshold exists.
  void addSelection(int observableIndex,
                    std::string const& selectionName,
                    int mode)
  {
    switch (mode) {
      case -1: // cut is optional and we store a bit for it
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{1}, limits::LimitType::kEqual, false, false, true);
        break;
      case 0: // cut is disabled; initialize with empty vector so evaluation bails out early
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{}, limits::LimitType::kEqual, false, false, false);
        break;
      case 1: // mandatory cut; only one threshold so the most permissive bit is skipped and no extra bit is stored
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{1}, limits::LimitType::kEqual, true, true, false);
        break;
      case 2: // pass through mode; cut is neither minimal nor optional and we store all bits
        mSelectionContainers.at(observableIndex) = SelectionContainer<T, BitmaskType>(selectionName, std::vector<T>{1}, limits::LimitType::kEqual, false, false, false);
        break;
      default:
        LOG(fatal) << "Invalid switch for boolean selection";
    }
    init(observableIndex);
  }

  /// \brief Update the limits of a function-based selection for a specific observable.
  /// \param observable Index of the observable.
  /// \param value Value at which to re-evaluate the selection functions.
  void updateLimits(int observable, T value)
  {
    mSelectionContainers.at(observable).updateLimits(value);
  }

  /// \brief Reset the internal bitmask and evaluation flags before processing a new candidate.
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

  /// \brief Reset the selection container for a single observable.
  /// \param observableIndex Index of the observable to reset.
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
    // if a previous observable already failed a minimal selection,
    // there is no point in evaluating further observables
    if (!mPassesMinimalSelections) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(value);
    // if any minimal selection is not fulfilled, the candidate is rejected
    if (mHasMinimalSelection) {
      if (!mSelectionContainers.at(observableIndex).passesAsMinimalCut()) {
        mPassesMinimalSelections = false;
      }
    }
    // if any optional selection is fulfilled, the candidate is accepted
    if (mHasOptionalSelection) {
      if (mSelectionContainers.at(observableIndex).passesAsOptionalCut()) {
        mPassesOptionalSelections = true;
      }
    }
  }

  /// \brief Evaluate a single observable against its configured selections.
  /// \param observableIndex Index of the observable.
  /// \param values Vector of values of the observable.
  void evaluateObservable(int observableIndex, std::vector<T> values)
  {
    // if there are no selections configured, bail out
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }
    // if a previous observable already failed a minimal selection,
    // there is no point in evaluating further observables
    if (!mPassesMinimalSelections) {
      return;
    }
    // set bitmask for given observable
    mSelectionContainers.at(observableIndex).evaluate(values);
    // if any minimal selection is not fulfilled, the candidate is rejected
    if (mHasMinimalSelection) {
      if (!mSelectionContainers.at(observableIndex).passesAsMinimalCut()) {
        mPassesMinimalSelections = false;
      }
    }
    // if any optional selection is fulfilled, the candidate is accepted
    if (mHasOptionalSelection) {
      if (mSelectionContainers.at(observableIndex).passesAsOptionalCut()) {
        mPassesOptionalSelections = true;
      }
    }
  }

  /// \brief Add comments to the selections of a specific observable.
  /// \param observableIndex Index of the observable.
  /// \param comments Vector of comment strings, one per selection threshold.
  void addComments(int observableIndex, std::vector<std::string> const& comments) { mSelectionContainers.at(observableIndex).addComments(comments); }

  /// \brief Check whether all required and optional cuts are passed.
  /// \return True if all minimal cuts pass and, if optional cuts are present, at least one of them passes.
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
    // if there are no minimal or optional selections, we pass let it pass with true
    return true;
  }

  /// \brief Check whether the optional selection for a specific observable is passed.
  /// \param observableIndex Index of the observable.
  /// \return True if at least one optional selection for this observable is fulfilled.
  bool passesOptionalSelection(int observableIndex) const
  {
    return mSelectionContainers.at(observableIndex).passesAsOptionalCut();
  }

  /// \brief Assemble the global selection bitmask from all individual observable selections.
  /// \tparam HistName Name of the histogram used to track selection statistics.
  template <const char* HistName>
  void assembleBitmask()
  {
    mHistRegistry->fill(HIST(HistName), mNSelection);
    // if the required selections are not passed, clear the bitmask and return early
    if (!this->passesAllRequiredSelections()) {
      mFinalBitmask.reset();
      return;
    }
    mHistRegistry->fill(HIST(HistName), mNSelection + 1);

    int binCenter = 0;
    // assemble the final bitmask by shifting each container's bitmask to its offset and OR-ing it in
    for (auto const& selectionContainer : mSelectionContainers) {
      // skip observables with no configured selections
      if (selectionContainer.isEmpty()) {
        continue;
      }
      // shift the container's bitmask to its offset and merge
      mFinalBitmask |= (selectionContainer.getBitmask() << selectionContainer.getOffset());

      for (std::size_t j = 0; j < selectionContainer.getNSelections(); ++j) {
        if (j == 0 && selectionContainer.skipMostPermissiveBit()) {
          // if the most permissive bit is skipped, this means this cut always applies
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
  /// \return The combined selection bitmask for all observables.
  BitmaskType getBitmask() const { return static_cast<BitmaskType>(mFinalBitmask.to_ullong()); }

  /// \brief Retrieve the bitmask for a single observable as an integer value.
  /// \param observableIndex Index of the observable.
  /// \return The selection bitmask for the specified observable.
  BitmaskType getBitmask(int observableIndex) const { return static_cast<BitmaskType>(mSelectionContainers.at(observableIndex).getBitmask().to_ullong()); }

  /// \brief Manually set the bitmask for a specific observable.
  /// \tparam R Integer type of the bitmask value.
  /// \param observableIndex Index of the observable.
  /// \param bitmask Bitmask value to set.
  template <typename R>
  void setBitmask(int observableIndex, R bitmask)
  {
    mSelectionContainers.at(observableIndex).setBitmask(bitmask);
  }

  /// \brief Retrieve the loosest (most permissive) selection threshold for a specific observable.
  /// \param observableIndex Index of the observable.
  /// \return The loosest threshold value configured for this observable.
  T getLoosestSelection(int observableIndex) const { return mSelectionContainers.at(observableIndex).getLoosestSelection(); }

  /// \brief Print the full configuration of all selections to the log.
  /// \param objectName Name of the object owning this selection (used as label in the log output).
  void printSelections(const std::string& objectName) const
  {
    LOG(info) << "Printing Configuration of " << objectName;
    for (size_t idx = 0; idx < mSelectionContainers.size(); ++idx) {
      const auto& container = mSelectionContainers.at(idx);
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

      for (std::size_t j = 0; j < container.getNSelections(); ++j) {

        std::stringstream line;
        std::string sel = container.getValueAsString(j);
        std::string comment = container.getComment(j);

        line << "    " << std::left << std::setw(25) << sel;

        if (j == 0 && container.skipMostPermissiveBit()) {
          line << "-> most permissive bit not saved";
        } else {
          int bit = container.getOffset() + (j - (container.skipMostPermissiveBit() ? 1 : 0));
          line << "-> Bit: 0x" << std::hex << std::uppercase << (1ULL << bit) << std::dec;
        }

        if (!comment.empty()) {
          line << " (" << comment << ")";
        }
        LOG(info) << line.str();
      }
      LOG(info) << "";
    }
    LOG(info) << "Number of occupied bits: " << mNSelectionBits << " / " << sizeof(BitmaskType) * CHAR_BIT;
    LOG(info) << "Printing done";
  }

  /// \brief Initialize histograms and set bitmask offsets for all configured observables.
  /// \tparam HistName Name of the histogram to create in the registry.
  /// \param registry Pointer to the histogram registry.
  template <const char* HistName>
  void setupContainers(o2::framework::HistogramRegistry* registry)
  {
    mHistRegistry = registry;
    // create histogram with one bin per selection, plus two summary bins (all analyzed, all passed)
    int nBins = mNSelection + 2;
    mHistRegistry->add(HistName, "; Selection Bits; Entries", o2::framework::HistType::kTH1F, {{nBins, -0.5, nBins - 0.5}});

    int binIndex = 0;
    int offset = 0;
    for (std::size_t idx = 0; idx < mSelectionContainers.size(); ++idx) {
      auto& container = mSelectionContainers[idx];
      if (container.isEmpty()) {
        continue;
      }
      container.setOffset(offset);
      offset += container.getShift();
      for (std::size_t j = 0; j < container.getNSelections(); j++) {
        std::string label = container.getBinLabel(j);
        mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(binIndex + 1, label.c_str());
        binIndex++;
      }
    }
    mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(mNSelection + 1, "All analyzed");
    mHistRegistry->get<TH1>(HIST(HistName))->GetXaxis()->SetBinLabel(mNSelection + 2, "All passed");
  }

 protected:
  void init(int observableIndex)
  {
    // check if any selections are configured
    if (mSelectionContainers.at(observableIndex).isEmpty()) {
      return;
    }

    // track the number of occupied bits and total selections
    mNSelectionBits += mSelectionContainers.at(observableIndex).getShift();
    mNSelection += mSelectionContainers.at(observableIndex).getNSelections();

    // check if any selection is minimal
    if (mSelectionContainers.at(observableIndex).isMinimalCut()) {
      mHasMinimalSelection = true;
    }
    // check if selection is optional
    if (mSelectionContainers.at(observableIndex).isOptionalCut()) {
      mHasOptionalSelection = true;
    }

    if (mNSelectionBits > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " number of bits are supported";
    }
  }

  o2::framework::HistogramRegistry* mHistRegistry = nullptr;
  std::array<SelectionContainer<T, BitmaskType>, NumObservables> mSelectionContainers = {}; ///< Array of selection containers, one per observable
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mFinalBitmask = {};                           ///< Assembled bitmask combining all observable selections
  std::size_t mNSelectionBits = 0;                                                          ///< Number of bits occupied in the bitmask (excludes skipped most-permissive bits)
  std::size_t mNSelection = 0;                                                              ///< Total number of configured selection thresholds across all observables
  bool mHasMinimalSelection = false;                                                        ///< True if at least one observable has a mandatory (minimal) cut configured
  bool mPassesMinimalSelections = true;                                                     ///< True if all mandatory (minimal) cuts have been passed so far
  bool mHasOptionalSelection = false;                                                       ///< True if at least one observable has an optional cut configured
  bool mPassesOptionalSelections = false;                                                   ///< True if at least one optional cut has been passed
};
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_BASESELECTION_H_
