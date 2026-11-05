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

#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>

#include <TF1.h>

#include <algorithm>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::analysis::femto
{

/// Limit type for selections
namespace limits
{
enum LimitType { kUpperLimit,            ///< simple upper limit, e.g. p_T < 1 GeV/c
                 kAbsUpperLimit,         ///< upper limit on the absolute value, e.g. |eta| < 0.8
                 kLowerLimit,            ///< simple lower limit, e.g. p_T > 0.2 GeV/c
                 kAbsLowerLimit,         ///< lower limit on the absolute value, e.g. |DCA_xy| > 0.05 cm
                 kUpperFunctionLimit,    ///< upper limit given by a function, e.g. DCA_xy < f(pt)
                 kAbsUpperFunctionLimit, ///< upper limit on the absolute value given by a function, e.g. |DCA_xy| < f(pt)
                 kLowerFunctionLimit,    ///< lower limit given by a function, e.g. DCA_xy > f(pt)
                 kAbsLowerFunctionLimit, ///< lower limit on the absolute value given by a function, e.g. |DCA_xy| > f(pt)
                 kEqual,                 ///< value must equal a fixed threshold, e.g. sign == 1
                 kEqualArray,            ///< each element of a value array must equal the corresponding threshold
                 kInsideRange,           ///< value must fall within [lower, upper]; either bound can be omitted (open interval)
                 kOutsideRange,          ///< value must fall outside [lower, upper]; fully-open range (both bounds omitted) disables the cut
                 kLimitTypeLast
};

inline const std::unordered_map<LimitType, std::string> limitTypeAsStrings = {
  {kUpperLimit, "Upper Limit"},
  {kAbsUpperLimit, "Absolute Upper Limit"},
  {kLowerLimit, "Lower Limit"},
  {kAbsLowerLimit, "Absolute Lower Limit"},
  {kUpperFunctionLimit, "Upper Function Limit"},
  {kAbsUpperFunctionLimit, "Absolute Upper Function Limit"},
  {kLowerFunctionLimit, "Lower Function Limit"},
  {kAbsLowerFunctionLimit, "Absolute Lower Function Limit"},
  {kEqual, "Equal"},
  {kEqualArray, "EqualArray"},
  {kInsideRange, "Inside Range"},
  {kOutsideRange, "Outside Range"},
  {kLimitTypeLast, "Last Limit Type"},
};

}; // namespace limits

namespace selectioncontainer
{

// strings for selection/filter histogram
constexpr std::string_view SectionDelimiter = ":::";
constexpr std::string_view RangeDelimiter = ";";
constexpr std::string_view ValueDelimiter = "___";
constexpr std::string_view NoValue = "X";

/// \class SelectionContainer
/// \brief Stores and evaluates multiple selection thresholds for a single observable.
///
/// Selections can be based on static threshold values or dynamically evaluated TF1 functions.
/// Thresholds are sorted from most permissive to most restrictive so that evaluation can bail
/// out early once a threshold fails.
///
/// \tparam T Data type for selection values (mostly floats).
/// \tparam BitmaskType Integer type used for bitmask storage (e.g., uint8_t, uint32_t).
template <typename T, typename BitmaskType>
class SelectionContainer
{
 public:
  /// \brief Default constructor.
  SelectionContainer() = default;

  /// \brief Destructor.
  ~SelectionContainer() = default;

  /// \brief Constructor for static value-based selections.
  /// \param selectionName Name of the observable this container manages.
  /// \param selectionValues Vector of threshold values.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the most permissive threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this selection is mandatory (candidate is rejected if it fails).
  /// \param isOptionalCut Whether this selection is optional (candidate is accepted if any optional cut passes).
  SelectionContainer(std::string const& selectionName,
                     std::vector<T> const& selectionValues,
                     limits::LimitType limitType,
                     bool skipMostPermissiveBit,
                     bool isMinimalCut,
                     bool isOptionalCut)
    : mSelectionName(selectionName),
      mSelectionValues(selectionValues),
      mLimitType(limitType),
      mSkipMostPermissiveBit(skipMostPermissiveBit),
      mIsMinimalCut(isMinimalCut),
      mIsOptionalCut(isOptionalCut)
  {
    if (mSelectionValues.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Selection '" << mSelectionName << "': too many selections for a single observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    if (isMinimalCut && isOptionalCut) {
      LOG(fatal) << "Selection '" << mSelectionName << "': a selection cannot be both minimal and optional";
    }
    // values for selection are not necessarily ordered correctly
    sortSelections();
  }

  /// \brief Constructor for range-based selections defined as "lower;upper" strings.
  ///        Either bound may be omitted to represent an open interval, e.g. ";1.0" means value < 1.0.
  ///        A fully-open range (both bounds omitted, i.e. ";") is exempt: for kInsideRange it means
  ///        "no constraint" (every value is inside); for kOutsideRange it means "no rejection" (every
  ///        value counts as outside) -- the natural reading breaks down for kOutsideRange otherwise,
  ///        since an infinitely wide excluded region would reject everything instead of nothing.
  /// \param selectionName Name of the observable this container manages.
  /// \param rangeStrings Vector of range strings, each of the form "lower;upper".
  /// \param skipMostPermissiveBit Whether to skip the most permissive threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this selection is mandatory (candidate is rejected if it fails).
  /// \param isOptionalCut Whether this selection is optional (candidate is accepted if any optional cut passes).
  /// \param limitType Either kInsideRange (value must fall inside the interval) or kOutsideRange
  ///                  (value must fall outside it, e.g. rejecting a competing mass hypothesis).
  SelectionContainer(std::string const& selectionName,
                     std::vector<std::string> const& rangeStrings,
                     bool skipMostPermissiveBit,
                     bool isMinimalCut,
                     bool isOptionalCut,
                     limits::LimitType limitType = limits::kInsideRange)
    : mSelectionName(selectionName),
      mLimitType(limitType),
      mSkipMostPermissiveBit(skipMostPermissiveBit),
      mIsMinimalCut(isMinimalCut),
      mIsOptionalCut(isOptionalCut)
  {
    if (mLimitType != limits::kInsideRange && mLimitType != limits::kOutsideRange) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range-string selections only support kInsideRange or kOutsideRange";
    }
    if (rangeStrings.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Selection '" << mSelectionName << "': too many selections for a single observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    if (isMinimalCut && isOptionalCut) {
      LOG(fatal) << "Selection '" << mSelectionName << "': a selection cannot be both minimal and optional";
    }
    for (auto const& rangeStr : rangeStrings) {
      T lower, upper;
      parseRangeString(rangeStr, lower, upper);
      mSelectionRanges.emplace_back(lower, upper);
    }
    // for kInsideRange, ranges are sorted from widest (most permissive) to tightest; for
    // kOutsideRange it's the opposite -- the narrowest excluded region is the most permissive
    // (fewer values get rejected)
    // in principle the ranges do not have to include each other, exepct this is configured as minimal cut, check is added at the end
    // to cover both cases, this also means we always check all ranges, when assembling the bit mask
    sortSelections();
    // init mSelectionValues to be the widths of the intervals
    for (std::size_t i = 0; i < mSelectionRanges.size(); i++) {
      mSelectionValues.push_back(mSelectionRanges[i].second - mSelectionRanges[i].first);
    }

    // for a minimal range cut, ranges must be strictly nested (each range contains all ranges that
    // sort after it) -- this is required for the single-bit shortcut in passesAsMinimalCut() to be
    // correct. For kInsideRange (sorted widest-first) that means range[i] must contain range[i+1];
    // for kOutsideRange (sorted narrowest-first) it's the reverse: range[i] must be contained in
    // range[i+1]. A fully-open ("reject nothing") entry is the empty interval [max, -max], which is
    // trivially contained in (or containing, for kInsideRange) anything, so this check already holds
    // for it without any special-casing.
    if (isMinimalCut) {
      for (std::size_t i = 0; i + 1 < mSelectionRanges.size(); i++) {
        bool const nested = (mLimitType == limits::kOutsideRange)
                              ? (mSelectionRanges[i].first >= mSelectionRanges[i + 1].first &&
                                 mSelectionRanges[i].second <= mSelectionRanges[i + 1].second)
                              : (mSelectionRanges[i].first <= mSelectionRanges[i + 1].first &&
                                 mSelectionRanges[i].second >= mSelectionRanges[i + 1].second);
        if (!nested) {
          LOG(fatal) << "Selection '" << selectionName << "': ranges for minimal cut"
                     << " are not nested. Range [" << mSelectionRanges[i].first << ";" << mSelectionRanges[i].second
                     << "] is not " << (mLimitType == limits::kOutsideRange ? "contained in" : "containing")
                     << " [" << mSelectionRanges[i + 1].first << ";" << mSelectionRanges[i + 1].second << "]";
        }
      }
    }
  }

  /// \brief Constructor for function-based dynamic selections.
  /// \param selectionName Name of the observable this container manages.
  /// \param lowerLimit Lower bound of the TF1 domain.
  /// \param upperLimit Upper bound of the TF1 domain.
  /// \param functions Vector of strings defining TF1 threshold functions.
  /// \param limitType Type of limit (from limits::LimitType).
  /// \param skipMostPermissiveBit Whether to skip the most permissive threshold when assembling the bitmask.
  /// \param isMinimalCut Whether this selection is mandatory (candidate is rejected if it fails).
  /// \param isOptionalCut Whether this selection is optional (candidate is accepted if any optional cut passes).
  SelectionContainer(std::string const& selectionName,
                     T lowerLimit,
                     T upperLimit,
                     std::vector<std::string> const& functions,
                     limits::LimitType limitType,
                     bool skipMostPermissiveBit,
                     bool isMinimalCut,
                     bool isOptionalCut)
    : mSelectionName(selectionName),
      mLimitType(limitType),
      mSkipMostPermissiveBit(skipMostPermissiveBit),
      mIsMinimalCut(isMinimalCut),
      mIsOptionalCut(isOptionalCut)
  {
    if (functions.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Selection '" << mSelectionName << "': too many selections for a single observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    if (isMinimalCut && isOptionalCut) {
      LOG(fatal) << "Selection '" << mSelectionName << "': a selection cannot be both minimal and optional";
    }
    for (std::size_t i = 0; i < functions.size(); i++) {
      mSelectionFunctions.emplace_back((mSelectionName + std::to_string(i)).c_str(), functions.at(i).c_str(), lowerLimit, upperLimit);
    }
    // functions are not necessarily ordered correctly;
    // use the midpoint of the domain to establish their order.
    // this relies on the user ensuring the ordering is consistent across the whole interval.
    T midPoint = (lowerLimit + upperLimit) / 2.;
    sortFunctions(midPoint);
    // initialize threshold values to the functions evaluated at the midpoint
    for (std::size_t i = 0; i < functions.size(); i++) {
      mSelectionValues.push_back(mSelectionFunctions.at(i).Eval(midPoint));
    }
  }

  /// \brief Attach comments to the selection thresholds, one per threshold in the same order as the input values.
  /// \param comments Vector of comment strings.
  void addComments(std::vector<std::string> const& comments)
  {
    // note: threshold values may be reordered by sortSelections() or sortFunctions(),
    // so comments must be provided in the already-sorted order
    if (comments.size() != getNSelections()) {
      LOG(fatal) << "Selection '" << mSelectionName << "': number of comments and number of selections are inconsistent";
    }
    mComments = comments;
  }

  /// \brief Attach the same comment to every threshold of this selection.
  /// \param comment Comment string applied to all thresholds.
  void addComments(std::string const& comment)
  {
    mComments = std::vector<std::string>(getNSelections(), comment);
  }

  /// \brief Get comments attached to the selection thresholds.
  /// \return Vector of comment strings.
  [[nodiscard]] std::string getComment(int selectionIndex) const
  {
    if (mComments.empty()) {
      return std::string{""};
    }
    return mComments.at(selectionIndex);
  }

  /// \brief Get the name of this selection.
  /// \return Selection name string.
  [[nodiscard]] std::string const& getSelectionName() const { return mSelectionName; }

  /// \brief Update threshold values by re-evaluating the internal TF1 functions at a given point.
  /// \param value Input value at which to evaluate the functions.
  void updateLimits(T value)
  {
    // functions are already sorted, so evaluate in the same order as mSelectionValues
    for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
      mSelectionValues.at(i) = mSelectionFunctions.at(i).Eval(value);
    }
  }

  /// \brief Evaluate which selection thresholds are passed for a given observable value.
  /// \param value Value of the observable to evaluate.
  void evaluate(T value)
  {
    // reset the bitmask before evaluating
    mBitmask.reset();

    // switch on limit type once outside the loop;
    // thresholds are sorted from most permissive to most restrictive,
    // so we can break early as soon as one comparison fails
    switch (mLimitType) {
      case (limits::kUpperLimit):
      case (limits::kUpperFunctionLimit):
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (value <= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            break;
          }
        }
        break;
      case (limits::kAbsUpperLimit):
      case (limits::kAbsUpperFunctionLimit):
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (std::abs(value) <= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            break;
          }
        }
        break;
      case (limits::kLowerLimit):
      case (limits::kLowerFunctionLimit):
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (value >= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            break;
          }
        }
        break;
      case (limits::kAbsLowerLimit):
      case (limits::kAbsLowerFunctionLimit):
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (std::abs(value) >= mSelectionValues.at(i)) {
            mBitmask.set(i);
          } else {
            break;
          }
        }
        break;
      case (limits::kInsideRange):
        // ranges are sorted widest-first but a narrower range passing does not imply a wider one fails (no check on boundries),
        // so all ranges must be checked explicitly — no early exit
        for (std::size_t i = 0; i < mSelectionRanges.size(); i++) {
          if (value >= mSelectionRanges.at(i).first && value <= mSelectionRanges.at(i).second) {
            mBitmask.set(i);
          }
        }
        break;
      case (limits::kOutsideRange):
        // ranges are sorted narrowest-first (most permissive -- smallest excluded region); same
        // no-early-exit reasoning as kInsideRange applies, just with the pass condition inverted.
        // A fully-open range string (";") is parsed as the empty interval [max, -max] (see
        // parseRangeString()), so this plain check already evaluates true for every value -- no
        // special-casing needed here to express "reject nothing".
        for (std::size_t i = 0; i < mSelectionRanges.size(); i++) {
          if (value < mSelectionRanges.at(i).first || value > mSelectionRanges.at(i).second) {
            mBitmask.set(i);
          }
        }
        break;
      case (limits::kEqual):
        // kEqual has no natural ordering, so all thresholds must be checked and we cannot bail early
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (std::fabs(value - mSelectionValues.at(i)) < constants::math::Epsilon) {
            mBitmask.set(i);
          }
        }
        break;
      default:
        LOG(warn) << "Limit type not known, no selection is applied";
    }
  }

  /// \brief Evaluate which selection thresholds are passed for a vector of observable values.
  ///        Only kEqualArray is supported for now; each element is compared against the corresponding threshold.
  /// \param values Values of the observable to evaluate.
  void evaluate(std::vector<T> const& values)
  {
    if (values.size() != mSelectionValues.size()) {
      LOG(fatal) << "Selection '" << mSelectionName << "': wrong number of values have been passed";
    }
    // reset the bitmask before evaluating
    mBitmask.reset();

    switch (mLimitType) {
      case (limits::kEqualArray):
        for (std::size_t i = 0; i < mSelectionValues.size(); i++) {
          if (std::fabs(values.at(i) - mSelectionValues.at(i)) < constants::math::Epsilon) {
            mBitmask.set(i);
          }
        }
        break;
      default:
        LOG(warn) << "Limit type not known, no selection is applied";
    }
  }

  /// \brief Retrieve the bitmask indicating which thresholds were passed.
  ///        If mSkipMostPermissiveBit is set, the bit for the loosest threshold is removed by shifting.
  /// \return Bitset representing passed selections.
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> getBitmask() const
  {
    if (!mSkipMostPermissiveBit) {
      return mBitmask;
    }
    // remove the first (most permissive) bit since it corresponds to the minimal selection and is always true
    return mBitmask >> 1;
  }

  /// \brief Manually set the internal bitmask.
  /// \tparam R Integer type of the bitmask value.
  /// \param bitmask Bitmask value to set.
  template <typename R>
  void setBitmask(R bitmask)
  {
    mBitmask = std::bitset<sizeof(BitmaskType) * CHAR_BIT>(bitmask);
  }

  /// \brief Retrieve the raw (unshifted) internal bitmask, including the most-permissive-threshold
  ///        bit even if mSkipMostPermissiveBit is set. Used when two separate evaluate() calls for
  ///        the same observable need to be combined externally (e.g. logical AND) before the
  ///        result is written back via setBitmask() -- unlike getBitmask(), this never shifts, so
  ///        it stays the correct counterpart to setBitmask() (which also never shifts).
  /// \return Bitset exactly as currently stored.
  [[nodiscard]] std::bitset<sizeof(BitmaskType) * CHAR_BIT> const& getRawBitmask() const { return mBitmask; }

  /// \brief Reset the internal bitmask to zero.
  void reset() { mBitmask.reset(); }

  /// \brief Check whether the mandatory (minimal) cut condition is fulfilled.
  /// \return True if the minimal selection passes or if this container is not marked as a minimal cut.
  [[nodiscard]] bool passesAsMinimalCut() const
  {
    if (mIsMinimalCut) {
      // if any bit is set the loosest threshold passed; since thresholds are ordered,
      // the loosest passing implies all looser ones also pass
      return mBitmask.test(0);
    }
    // not a minimal cut — return true by default so it does not block the candidate
    return true;
  }

  /// \brief Check whether any optional cut is fulfilled.
  /// \return True if at least one optional threshold is passed, false if this container is not marked as optional.
  [[nodiscard]] bool passesAsOptionalCut() const
  {
    if (mIsOptionalCut) {
      // if any bit is set the loosest threshold passed
      return mBitmask.any();
    }
    // not an optional cut — return false by default so it does not accidentally accept the candidate
    return false;
  }

  /// \brief Get the loosest (most permissive) selection threshold value.
  /// \return The first element of the sorted threshold vector.
  T getLoosestSelection() const
  {
    if (mSelectionValues.empty()) {
      LOG(fatal) << "Selection '" << mSelectionName << "': no selections configured";
    }
    return mSelectionValues.at(0);
  }

  /// \brief Check whether any selection thresholds are configured.
  ///        For function-based selections, mSelectionValues is always populated (initialised at the midpoint),
  ///        so this check is safe for both static and function-based containers.
  /// \return True if no thresholds are configured.
  [[nodiscard]] bool isEmpty() const { return mSelectionValues.empty(); }

  /// \brief Get the number of bits this container contributes to the global bitmask.
  ///        If the most permissive bit is skipped, the contribution is reduced by one.
  /// \return Number of bits to add to the global bitmask offset.
  [[nodiscard]] int getShift() const
  {
    if (mSelectionValues.empty()) {
      return 0;
    }
    if (mSkipMostPermissiveBit) {
      return static_cast<int>(mSelectionValues.size() - 1);
    }
    return static_cast<int>(mSelectionValues.size());
  }

  /// \brief Set the bit offset of this container within the global bitmask.
  /// \param offset Bit position at which this container's bits start.
  void setOffset(int offset) { mOffset = offset; }

  /// \brief Get the bit offset of this container within the global bitmask.
  /// \return Bit offset.
  [[nodiscard]] int getOffset() const { return mOffset; }

  /// \brief Get the total number of configured selection thresholds.
  /// \return Number of thresholds.
  [[nodiscard]] std::size_t getNSelections() const { return mSelectionValues.size(); }

  /// \brief Build a histogram bin label string encoding the full configuration of a single threshold.
  /// \param selectionIndex Index of the threshold within this container.
  /// \return Encoded label string.
  /// \brief Describes an open (sentinel, i.e. omitted-by-the-user) side of a range, if any.
  /// For kOutsideRange, a single omitted side flips which half of the "value < lower ||
  /// value > upper" OR-condition stays active, so the resulting satisfied/reject direction
  /// is spelled out explicitly rather than named generically -- omitting the upper bound
  /// there means the lower-bound inequality alone remains, i.e. "reject value < lower".
  /// \param range Lower/upper bound pair for a kInsideRange/kOutsideRange threshold.
  /// \return A human-readable description of the effective cut, or "" if both sides are bounded.
  [[nodiscard]] std::string getOpenSideComment(std::pair<T, T> const& range) const
  {
    bool const lowerOpen = std::abs(range.first) == std::numeric_limits<T>::max();
    bool const upperOpen = std::abs(range.second) == std::numeric_limits<T>::max();
    if (!lowerOpen && !upperOpen) {
      return "";
    }
    if (lowerOpen && upperOpen) {
      return "no cut";
    }
    std::ostringstream oss;
    if (mLimitType == limits::kInsideRange) {
      oss << "accept value " << (upperOpen ? ">= " : "<= ") << (upperOpen ? range.first : range.second);
    } else {
      oss << "reject value " << (upperOpen ? "< " : "> ") << (upperOpen ? range.first : range.second);
    }
    return oss.str();
  }

  [[nodiscard]] std::string getBinLabel(int selectionIndex) const
  {
    std::ostringstream oss;

    // Determine value string
    std::string valueStr;
    std::string autoComment;

    if (mLimitType == limits::kInsideRange || mLimitType == limits::kOutsideRange) {
      // Print the actual lower;upper interval -- kept as real, parseable numbers (not a word like
      // "no cut") since this Value field is read by cutculator.py; an open side (lower and/or
      // upper bound left at the sentinel, i.e. omitted by the user) is instead flagged via the
      // Comment field below.
      const auto& range = mSelectionRanges.at(selectionIndex);
      std::ostringstream rangeStream;
      rangeStream << range.first << RangeDelimiter << range.second;
      valueStr = rangeStream.str();
      autoComment = getOpenSideComment(range);
    } else if (mSelectionFunctions.empty()) {
      valueStr = std::to_string(mSelectionValues.at(selectionIndex));
    } else {
      valueStr = mSelectionFunctions.at(selectionIndex).GetExpFormula().Data();
    }

    // combine rather than override -- e.g. TOF PID selections attach a "minMomTofX = ..." comment
    // to every threshold via addComments(), independent of whether a given threshold's range
    // also happens to be open-sided
    std::string const userComment = mComments.empty() ? "" : mComments.at(selectionIndex);
    std::string comment = userComment;
    if (!autoComment.empty()) {
      comment = userComment.empty() ? autoComment : autoComment + "; " + userComment;
    }

    oss << "SelectionName" << ValueDelimiter << mSelectionName << SectionDelimiter
        << "LimitType" << ValueDelimiter << getLimitTypeAsString() << SectionDelimiter
        << "MinimalCut" << ValueDelimiter << (mIsMinimalCut ? "1" : "0") << SectionDelimiter
        << "SkipMostPermissiveBit" << ValueDelimiter << (mSkipMostPermissiveBit ? "1" : "0") << SectionDelimiter
        << "OptionalCut" << ValueDelimiter << (mIsOptionalCut ? "1" : "0") << SectionDelimiter
        << "Shift" << ValueDelimiter << getShift() << SectionDelimiter
        << "Offset" << ValueDelimiter << mOffset << SectionDelimiter
        << "Value" << ValueDelimiter << valueStr << SectionDelimiter
        << "BitPosition" << ValueDelimiter << (mSkipMostPermissiveBit ? (selectionIndex == 0 ? NoValue : std::to_string(mOffset + selectionIndex - 1)) : std::to_string(mOffset + selectionIndex)) << SectionDelimiter
        << "Comment" << ValueDelimiter << (comment.empty() ? NoValue : comment);
    return oss.str();
  }

  [[nodiscard]] std::string getValueAsString(int selectionIndex) const
  {
    if (this->isEmpty()) {
      return std::string{"No value configured"};
    }
    if (!mSelectionFunctions.empty()) {
      return std::string(mSelectionFunctions.at(selectionIndex).GetExpFormula().Data());
    }
    if (!mSelectionRanges.empty()) {
      const auto& range = mSelectionRanges.at(selectionIndex);
      bool const lowerOpen = std::abs(range.first) == std::numeric_limits<T>::max();
      bool const upperOpen = std::abs(range.second) == std::numeric_limits<T>::max();
      if (lowerOpen && upperOpen) {
        return std::string{"no cut"};
      }
      std::ostringstream oss;
      oss << (lowerOpen ? "-inf" : std::to_string(range.first).c_str())
          << ";"
          << (upperOpen ? "+inf" : std::to_string(range.second).c_str());
      return oss.str();
    }
    return std::to_string(mSelectionValues.at(selectionIndex));
  }

  /// \brief Get the global bit position of a threshold within the final bitmask.
  ///        Calling this for the skipped most-permissive threshold is a fatal error.
  /// \param selectionIndex Index of the threshold within this container.
  /// \return Global bit position.
  [[nodiscard]] int getBitPosition(int selectionIndex) const
  {
    if (selectionIndex == 0 && mSkipMostPermissiveBit) {
      LOG(fatal) << "Selection '" << mSelectionName << "': trying to access the bit position of a skipped selection. Breaking...";
      return -1;
    }
    if (mSkipMostPermissiveBit) {
      return mOffset + selectionIndex - 1;
    }
    return mOffset + selectionIndex;
  }

  /// \brief Get the string representation of the configured limit type.
  /// \return Human-readable limit type name.
  [[nodiscard]] std::string getLimitTypeAsString() const { return limits::limitTypeAsStrings.at(mLimitType); }

  /// \brief Get the configured static threshold values.
  /// \return Const reference to the vector of threshold values.
  std::vector<T> const& getSelectionValues() const { return mSelectionValues; }

  /// \brief Get the configured TF1 threshold functions.
  /// \return Const reference to the vector of TF1 functions.
  [[nodiscard]] std::vector<TF1> const& getSelectionFunction() const { return mSelectionFunctions; }

  /// \brief Check whether this container is marked as a mandatory (minimal) cut.
  /// \return True if this is a minimal cut.
  [[nodiscard]] bool isMinimalCut() const { return mIsMinimalCut; }

  /// \brief Check whether this container is marked as an optional cut.
  /// \return True if this is an optional cut.
  [[nodiscard]] bool isOptionalCut() const { return mIsOptionalCut; }

  /// \brief Check whether the most permissive threshold bit is skipped when assembling the bitmask.
  /// \return True if the most permissive bit is skipped.
  [[nodiscard]] [[nodiscard]] bool skipMostPermissiveBit() const { return mSkipMostPermissiveBit; }

 private:
  /// \brief Sort static threshold values from most permissive to most restrictive based on the limit type.
  void sortSelections()
  {
    switch (mLimitType) {
      case (limits::kUpperLimit):
      case (limits::kAbsUpperLimit):
        std::sort(mSelectionValues.begin(), mSelectionValues.end(), [](T a, T b) { return a > b; });
        break;
      case (limits::kLowerLimit):
      case (limits::kAbsLowerLimit):
        std::sort(mSelectionValues.begin(), mSelectionValues.end(), [](T a, T b) { return a < b; });
        break;
      case (limits::kInsideRange):
        // sort by range width descending so the most permissive (widest) range comes first
        std::sort(mSelectionRanges.begin(), mSelectionRanges.end(),
                  [](std::pair<T, T> const& a, std::pair<T, T> const& b) {
                    return (a.second - a.first) > (b.second - b.first);
                  });
        break;
      case (limits::kOutsideRange):
        // opposite of kInsideRange: the most permissive (most likely to be satisfied) exclusion is
        // the narrowest one, since it rejects the fewest values. A fully-open ("reject nothing")
        // range parses to the empty interval [max, -max] (see parseRangeString()), whose "width" is
        // the most negative possible value, so it naturally sorts first without any special-casing.
        std::sort(mSelectionRanges.begin(), mSelectionRanges.end(),
                  [](std::pair<T, T> const& a, std::pair<T, T> const& b) {
                    return (a.second - a.first) < (b.second - b.first);
                  });
        break;
      default:
        break;
    }
  }

  /// \brief Sort selection functions from most permissive to most restrictive, evaluated at a given point.
  /// \param value Point at which to evaluate the functions for ordering.
  void sortFunctions(T value)
  {
    switch (mLimitType) {
      case (limits::kUpperFunctionLimit):
      case (limits::kAbsUpperFunctionLimit):
        std::sort(mSelectionFunctions.begin(), mSelectionFunctions.end(), [value](TF1 const& a, TF1 const& b) { return a.Eval(value) > b.Eval(value); });
        break;
      case (limits::kLowerFunctionLimit):
      case (limits::kAbsLowerFunctionLimit):
        std::sort(mSelectionFunctions.begin(), mSelectionFunctions.end(), [value](TF1 const& a, TF1 const& b) { return a.Eval(value) < b.Eval(value); });
        break;
      default:
        break;
    }
  }

  /// \brief Parse a range string of the form "lower;upper" into a lower and upper bound.
  ///        Either bound may be omitted (e.g. ";1.0" or "0.5;"), in which case
  ///        -/+ numeric_limits::max() is used respectively. For kOutsideRange, a fully-open string
  ///        (";", both bounds omitted) is the one exception: it is parsed as the empty interval
  ///        [max, -max] instead of [-max, max] -- every value then trivially fails the normal
  ///        "value >= lower && value <= upper" check, i.e. every value counts as outside, which is
  ///        what "reject nothing" needs. [-max, max] would instead mean every value is inside,
  ///        i.e. nothing is ever outside, which would reject everything -- the opposite intent.
  /// \param rangeStr Input string to parse.
  /// \param lower Output lower bound.
  /// \param upper Output upper bound.
  void parseRangeString(std::string const& rangeStr, T& lower, T& upper) const
  {
    auto pos = rangeStr.find(';');
    if (pos == std::string::npos) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range string '" << rangeStr << "' is missing ';' separator. Expected format: 'lower;upper'";
    }
    if (rangeStr.find(';', pos + 1) != std::string::npos) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range string '" << rangeStr << "' has more than one ';' separator. Expected format: 'lower;upper'";
    }
    std::string lowerStr = rangeStr.substr(0, pos);
    std::string upperStr = rangeStr.substr(pos + 1);

    if (lowerStr.empty() && upperStr.empty() && mLimitType == limits::kOutsideRange) {
      // kOutsideRange's "value < lower || value > upper" is the negation of kInsideRange's
      // "value >= lower && value <= upper" -- substituting the same [-max, +max] that makes the
      // latter trivially true for every value would make the former trivially FALSE for every
      // value instead (nothing is < -max or > +max), i.e. "reject everything" rather than the
      // intended "reject nothing". The empty interval [max, -max] is the encoding that keeps the
      // disjunction trivially true for every value instead, giving the correct "exempt" meaning.
      lower = std::numeric_limits<T>::max();
      upper = -std::numeric_limits<T>::max();
      return;
    }

    lower = lowerStr.empty() ? -std::numeric_limits<T>::max() : parseBound(rangeStr, lowerStr);
    upper = upperStr.empty() ? std::numeric_limits<T>::max() : parseBound(rangeStr, upperStr);

    if (lower >= upper) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range string '" << rangeStr << "' has lower bound >= upper bound";
    }
  }

  /// \brief Parse a single range bound and make sure the whole string was consumed, so that
  ///        malformed input such as "3;5" (extra ';') or "-3,5" (comma used as decimal
  ///        separator) is rejected instead of silently truncated by std::stod (which happily
  ///        parses just the numeric prefix "3" or "-3" and ignores the rest).
  /// \param rangeStr Full "lower;upper" string, only used for the error message.
  /// \param boundStr The individual bound substring (lower or upper) to parse.
  T parseBound(std::string const& rangeStr, std::string const& boundStr) const
  {
    std::size_t idx = 0;
    double value = 0.;
    try {
      value = std::stod(boundStr, &idx);
    } catch (std::exception const& e) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range string '" << rangeStr << "' contains bound '" << boundStr << "' that is not a valid number (" << e.what() << ")";
    }
    if (idx != boundStr.size()) {
      LOG(fatal) << "Selection '" << mSelectionName << "': range string '" << rangeStr << "' contains bound '" << boundStr
                 << "' with trailing characters ('" << boundStr.substr(idx) << "') that could not be parsed";
    }
    return static_cast<T>(value);
  }

  std::string mSelectionName;
  std::vector<T> mSelectionValues = {};                      ///< Threshold values, sorted from most permissive to most restrictive
  std::vector<TF1> mSelectionFunctions;                      ///< TF1 threshold functions (empty for static selections)
  std::vector<std::pair<T, T>> mSelectionRanges = {};        ///< Lower and upper bounds for kInsideRange/kOutsideRange selections, one pair per threshold
  limits::LimitType mLimitType = limits::kLimitTypeLast;     ///< Comparison type applied during evaluation
  bool mSkipMostPermissiveBit = false;                       ///< If true, the most permissive threshold does not occupy a bit in the global bitmask
  bool mIsMinimalCut = false;                                ///< If true, this selection is mandatory; failing it rejects the candidate
  bool mIsOptionalCut = false;                               ///< If true, this selection is optional; passing it accepts the candidate
  std::vector<std::string> mComments;                        ///< Optional comments per threshold, in the same order as mSelectionValues
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mBitmask = {}; ///< Bitmask indicating which thresholds were passed during the last evaluation
  int mOffset = 0;                                           ///< Bit offset of this container within the global bitmask
};
}; // namespace selectioncontainer

} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_SELECTIONCONTAINER_H_
