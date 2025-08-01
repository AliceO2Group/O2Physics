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

/// \file selectionContainer.h
/// \brief SelectionContainer - small container holding selections of an observable
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_
#define PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_

#include "CommonConstants/MathConstants.h"

#include "TF1.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <string>
#include <vector>

namespace o2::analysis::femtounited
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

std::unordered_map<LimitType, std::string> LimitTypeAsStrings = {
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

/// Simple class for storing selections of a single observable
/// \tparam T Data type used for the selection values (float/int/...)
/// \tparam BitmaskType Compute number of selections from BitmaskType (should be some unsigned integer)
template <typename T, typename BitmaskType>
class SelectionContainer
{
 public:
  /// Default constructor
  SelectionContainer() {}

  /// Constructor
  /// \param values Vector of values for the selection
  /// \param limitType Type of limit of the selection
  /// \param SkipLastBit Boolean whether to skip the last bit
  SelectionContainer(std::vector<T> const& SelectionValues, limits::LimitType limitType, bool SkipMostPermissiveBit, bool IsMinimalCut)
    : mSelectionValues(SelectionValues),
      mLimitType(limitType),
      mSkipMostPermissiveBit(SkipMostPermissiveBit),
      mIsMinimalCut(IsMinimalCut)
  {
    if (mSelectionValues.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    // for kEqual we can never skip the last bit
    if (limitType == limits::kEqual) {
      mSkipMostPermissiveBit = false;
    }
    // values for selection are not necessarily ordered correctly
    sortSelections();
  }

  /// Constructor
  /// \param baseName base name for TF1 object
  /// \param lowerLimit upper limit of the TF1 object
  /// \param upperLimit lower limit of the TF1 object
  /// \param function vector of strings for initializing of the TF1 object
  /// \param limitType Type of limit of the selection
  /// \param SkipLastBit Boolean whether to skip the last bit
  SelectionContainer(std::string const& baseName, T lowerLimit, T upperLimit, std::vector<std::string> const& functions, limits::LimitType limitType, bool skipMostPermissiveBit, bool IsMinimalCut)
    : mLimitType(limitType),
      mSkipMostPermissiveBit(skipMostPermissiveBit),
      mIsMinimalCut(IsMinimalCut)
  {
    if (functions.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    for (std::size_t i = 0; i < functions.size(); i++) {
      mSelectionFunctions.emplace_back((baseName + std::to_string(i)).c_str(), functions.at(i).c_str(), lowerLimit, upperLimit);
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

  /// Destructor
  virtual ~SelectionContainer() = default;

  /// Sort selections accroding to limit type
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

  // sort limit functions
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

  // update the selection limits depending on the passed function
  void updateLimits(T value)
  {
    // functions are ordered so just add the values in the same order
    for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
      mSelectionValues.at(i) = mSelectionFunctions.at(i).Eval(value);
    }
  }

  /// Check which selections are fulfilled
  /// \param observable Value of the variable to be checked
  void evaluate(T value)
  {
    // better safe than sorry and reset the bitmask before you evaluate and set minimal selection to true
    mBitmask.reset();
    // the values are ordered, from most loost to most tight, as soon as one comparison is not true, we can break out of the loop
    bool breakLoop = false;
    // iterate over all limits and set the corresponding bit if we pass the selection, otherwise break out as soon as we can
    // only break if the observable is used for the minimal selection, for example
    // an example, we configured |eta|<0.8 and  |nsigma_TPC_proton| < 3 and |nsigma_TPC_pion| < 3
    // all tracks need to be within |eta|<0.8, so we mark is as minimal selection
    // but the nsigma_TPC > 3 for proton hypothesis and we configured for protons and pions
    // the tracks would still be valid to use
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
          if (std::abs(value - mSelectionValues.at(i)) < constants::math::Epsilon) {
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

  /// Return the bitmask of the selections
  /// \return bitmask
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> getBitmask() const
  {
    // if we do not skip the last bit, return full bitmask
    // in the constructor we ensure that for kEqual we do not skip the most permissive bit
    if (mSkipMostPermissiveBit == false) {
      return mBitmask;
    } else {
      // for the other selections we can remove the first bit since it is the minimal selection and therefore always true
      return mBitmask >> 1;
    }
  }

  /// Check whether the minimal selection is fulfilled or not
  /// \return Whether the selection is fulfilled or not
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

  /// Return loosest selection value
  ///  The values are ordered, the loosest value will be the first one
  /// \return loosest selection
  T getLoosestSelection() const { return mSelectionValues.at(0); }

  /// Check whether  any selections are configured
  bool isEmpty() const { return mSelectionValues.empty(); }

  /// Get shift for final bitmask
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

  std::string getLimitTypeAsString() const { return limits::LimitTypeAsStrings[mLimitType]; }
  std::vector<T> getSelectionValues() const { return mSelectionValues; }
  bool isMinimalCut() const { return mIsMinimalCut; }
  bool skipMostPermissiveBit() const { return mSkipMostPermissiveBit; }

 private:
  std::vector<T> mSelectionValues = {};                      ///< Values used for the selection
  std::vector<TF1> mSelectionFunctions = {};                 ///< Function used for the selection
  limits::LimitType mLimitType;                              ///< Limit type of selection
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mBitmask = {}; ///< bitmask for the observable
  bool mSkipMostPermissiveBit = false;                       ///< whether to skip the last bit or not
  bool mIsMinimalCut = false;                                ///< whether to use this observable for minimal selection or not
};

} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_
