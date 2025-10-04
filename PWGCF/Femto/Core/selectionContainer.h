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

/// \file selectionContainer.h
/// \brief Defines the SelectionContainer class for managing selection criteria in analysis.
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_SELECTIONCONTAINER_H_
#define PWGCF_FEMTO_CORE_SELECTIONCONTAINER_H_

#include "CommonConstants/MathConstants.h"

#include "TF1.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{

/// Limit type for selections
namespace limits
{
enum LimitType { kUpperLimit,            ///< simple upper limit for the value, e.g. p_T < 1 GeV/c
                 kAbsUpperLimit,         ///< upper limit of the absolute value, e.g. |eta| < 0.8
                 kLowerLimit,            ///< simple lower limit for the value, e.g. p_T > 0.2 GeV/c
                 kAbsLowerLimit,         ///< lower limit of the absolute value, e.g. |DCA_xyz| > 0.05 cm
                 kEqual,                 ///< values need to be equal, e.g. sign = 1
                 kUpperFunctionLimit,    ///< simple upper limit of a function value, e.g. DCA_xy > f(pt)
                 kAbsUpperFunctionLimit, ///< upper limit of an absolute value given by a function, e.g. |DCA_xy| > f(pt)
                 kLowerFunctionLimit,    ///< simple lower limit of a function value, e.g. DCA_xy < f(pt)
                 kAbsLowerFunctionLimit  ///< lower limit of an absolute value given by a function, e.g. |DCA_xy| < f(pt)
};

std::unordered_map<LimitType, std::string> limitTypeAsStrings = {
  {kUpperLimit, "Upper Limit"},
  {kAbsUpperLimit, "Absolute Upper Limit"},
  {kLowerLimit, "Lower Limit"},
  {kAbsLowerLimit, "Absolute Lower Limit"},
  {kEqual, "Equal"},
  {kUpperFunctionLimit, "Upper Function Limit"},
  {kAbsUpperFunctionLimit, "Absolute Upper Function Limit"},
  {kLowerFunctionLimit, "Lower Function Limit"},
  {kAbsLowerFunctionLimit, "Absolute Lower Function Limit"}};

}; // namespace limits

/// \class SelectionContainer
/// \brief Class for storing and evaluating multiple selection thresholds for a single observable.
/// \tparam T Data type for selection values (mostly floats)
/// \tparam BitmaskType Type used for bitmask storage (e.g., uint8_t, uint32_t).
template <typename T, typename BitmaskType>
class SelectionContainer
{
 public:
  /// Default constructor
  SelectionContainer() {}

  /// \brief Constructor for static value-based selection.
  /// \param SelectionValues Vector of values for the selection.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param SkipMostPermissiveBit Whether to skip the most permissive bit in the bitmask.
  /// \param IsMinimalCut Whether this selection should be treated as a minimal required cut.
  SelectionContainer(std::vector<T> const& SelectionValues, limits::LimitType limitType, bool SkipMostPermissiveBit, bool IsMinimalCut)
    : mSelectionValues(SelectionValues),
      mLimitType(limitType),
      mSkipMostPermissiveBit(SkipMostPermissiveBit),
      mIsMinimalCut(IsMinimalCut)
  {
    if (mSelectionValues.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    // values for selection are not necessarily ordered correctly
    sortSelections();
  }

  /// \brief Constructor for function-based dynamic selection.
  /// \param baseName Base name for TF1 functions.
  /// \param lowerLimit Lower bound for TF1 domain.
  /// \param upperLimit Upper bound for TF1 domain.
  /// \param functions Vector of strings defining TF1 functions.
  /// \param limitType Type of limit.
  /// \param skipMostPermissiveBit Whether to skip the most permissive bit in the bitmask.
  /// \param IsMinimalCut Whether this selection should be treated as a minimal required cut.
  SelectionContainer(std::string const& baseName,
                     T lowerLimit,
                     T upperLimit,
                     std::vector<std::string> const& functions,
                     limits::LimitType limitType,
                     bool skipMostPermissiveBit,
                     bool IsMinimalCut)
    : mLimitType(limitType),
      mSkipMostPermissiveBit(skipMostPermissiveBit),
      mIsMinimalCut(IsMinimalCut)
  {
    if (functions.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    for (std::size_t i = 0; i < functions.size(); i++) {
      const std::string& func = functions.at(i);
      const std::string& safeFunc = func.empty() ? "0.1" : func; // in case string is empty, set to constant value of 0.1
      mSelectionFunctions.emplace_back((baseName + std::to_string(i)).c_str(), safeFunc.c_str(), lowerLimit, upperLimit);
    }
    // functions for selection are not necessarily ordered correctly
    // use value at midpoint to order them
    // here we rely on the user that the functions can be ordered like this over the whole interval
    T midPoint = (lowerLimit + upperLimit) / 2.;
    sortFunctions(midPoint);
    // initialize the values also to the midpoint
    for (std::size_t i = 0; i < functions.size(); i++) {
      mSelectionValues.push_back(mSelectionFunctions.at(i).Eval(midPoint));
    }
  }

  virtual ~SelectionContainer() = default;

  /// \brief Sort static selection values based on the limit type.
  void sortSelections()
  {
    switch (mLimitType) {
      case (limits::kUpperLimit):
      case (limits::kAbsUpperLimit):
        std::sort(mSelectionValues.begin(), mSelectionValues.end(), [](T a, T b) { return a > b; });
        break;
      case (limits::kLowerLimit):
      case (limits::kAbsLowerLimit):
      case (limits::kEqual):
        std::sort(mSelectionValues.begin(), mSelectionValues.end(), [](T a, T b) { return a < b; });
        break;
      default:
        break;
    }
  }

  /// \brief Sort selection functions based on evaluation at a given point.
  /// \param value Point at which to evaluate the functions for ordering.
  void sortFunctions(T value)
  {
    switch (mLimitType) {
      case (limits::kUpperFunctionLimit):
      case (limits::kAbsUpperFunctionLimit):
        std::sort(mSelectionFunctions.begin(), mSelectionFunctions.end(), [value](TF1 a, TF1 b) { return a.Eval(value) > b.Eval(value); });
        break;
      case (limits::kLowerFunctionLimit):
      case (limits::kAbsLowerFunctionLimit):
        std::sort(mSelectionFunctions.begin(), mSelectionFunctions.end(), [value](TF1 a, TF1 b) { return a.Eval(value) < b.Eval(value); });
        break;
      default:
        break;
    }
  }

  /// \brief Update selection limits using internal functions evaluated at a given value.
  /// \param value Input value to evaluate functions at.
  void updateLimits(T value)
  {
    // functions are ordered so just add the values in the same order
    for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
      mSelectionValues.at(i) = mSelectionFunctions.at(i).Eval(value);
    }
  }

  /// \brief Evaluate which selection criteria are fulfilled for a given value.
  /// \param value Value of the observable to evaluate.
  void evaluate(T value)
  {
    // better safe than sorry and reset the bitmask before you evaluate and set minimal selection to true
    mBitmask.reset();
    // the values are ordered, from most loose to most tight, as soon as one comparison is not true, we can break out of the loop
    bool breakLoop = false;
    // iterate over all limits and set the corresponding bit if we pass the selection, otherwise break out as soon as we can
    // only break if the observable is used for the minimal selection
    for (size_t i = 0; i < mSelectionValues.size(); i++) {
      switch (mLimitType) {
        case (limits::kUpperLimit):
        case (limits::kUpperFunctionLimit):
          if (value <= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kAbsUpperLimit):
        case (limits::kAbsUpperFunctionLimit):
          if (std::abs(value) <= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kLowerLimit):
        case (limits::kLowerFunctionLimit):
          if (value >= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kAbsLowerLimit):
        case (limits::kAbsLowerFunctionLimit):
          if (std::abs(value) >= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kEqual):
          // special case for kEqual since here we cannot really establish an order so we need to check all cases explicitly and we cannot bail early
          if (std::fabs(value - mSelectionValues.at(i)) < constants::math::Epsilon) {
            mBitmask.set(i);
          }
          break;
        default:
          breakLoop = true;
      }
      // bail early if a comparison fails
      // the values are ordered, so all following we also fail, there there is no point in contiuing
      if (breakLoop) {
        break;
      }
    }
  }

  /// \brief Retrieve the bitmask indicating which selections were passed.
  /// \return Bitset representing passed selections.
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> getBitmask() const
  {
    // if we do not skip the last bit, return full bitmask
    if (mSkipMostPermissiveBit == false) {
      return mBitmask;
    } else {
      // for the other selections we can remove the first bit since it is the minimal selection and therefore always true
      return mBitmask >> 1;
    }
  }
  template <typename R>
  void setBitmask(R bitmask)
  {
    mBitmask = std::bitset<sizeof(BitmaskType) * CHAR_BIT>(bitmask);
  }

  /// \brief Check whether the minimal cut condition is fulfilled.
  /// \return True if minimal selection is fulfilled, false otherwise.
  bool passesAsMinimalCut() const
  {
    if (mIsMinimalCut) {
      // check if loosest bit is set
      return mBitmask.test(0);
    } else {
      // if selection is not marked as a minimal cut, we return true by default
      return true;
    }
  }

  /// \brief Check whether any optional cuts are fulfilled.
  /// \return True if at least one optional cut is passed.
  bool passesAsOptionalCut() const
  {
    // if selection is marekd as minimal cut, we return false by default
    if (mIsMinimalCut) {
      return false;
    } else {
      // check if any bit is set
      return mBitmask.any();
    }
  }

  /// \brief Get the loosest (most permissive) selection value.
  /// \return First (loosest) selection value.
  T getLoosestSelection() const { return mSelectionValues.at(0); }

  /// \brief Check if there are any selection values configured.
  /// \return True if no selections are configured.
  bool isEmpty() const { return mSelectionValues.empty(); }

  /// \brief Get the number of bits to shift for the final bitmask.
  /// \return Number of bits to shift.
  int getShift() const
  {
    if (mSelectionValues.empty()) {
      return 0;
    }
    if (mSkipMostPermissiveBit) {
      return static_cast<int>(mSelectionValues.size() - 1);
    } else {
      return static_cast<int>(mSelectionValues.size());
    }
  }

  /// \brief Get string representation of the limit type.
  /// \return String name of the limit type.
  std::string getLimitTypeAsString() const { return limits::limitTypeAsStrings[mLimitType]; }

  /// \brief Get a copy of all selection values.
  /// \return Vector of selection values.
  std::vector<T> getSelectionValues() const { return mSelectionValues; }

  /// \brief Get a copy of all selection values.
  /// \return Vector of selection values.
  std::vector<TF1> getSelectionFunction() const { return mSelectionFunctions; }

  /// \brief Check if this container is marked as minimal cut.
  /// \return True if minimal cut, false otherwise.
  bool isMinimalCut() const { return mIsMinimalCut; }

  /// \brief Check whether the most permissive bit is skipped.
  /// \return True if skipped, false otherwise.
  bool skipMostPermissiveBit() const { return mSkipMostPermissiveBit; }

 private:
  std::vector<T> mSelectionValues = {};                      ///< Values used for the selection
  std::vector<TF1> mSelectionFunctions = {};                 ///< Function used for the selection
  limits::LimitType mLimitType;                              ///< Limit type of selection
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mBitmask = {}; ///< bitmask for the observable
  bool mSkipMostPermissiveBit = false;                       ///< whether to skip the last bit or not
  bool mIsMinimalCut = false;                                ///< whether to use this observable for minimal selection or not
};

} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_SELECTIONCONTAINER_H_
