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

/// \file FemtoUniverseSelection.h
/// \brief FemtoUniverseSelection - small generic class to do selections
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESELECTION_H_

namespace o2::analysis::femtoUniverse
{

namespace femtoUniverseSelection
{
/// Type of selection to be employed
enum SelectionType { kUpperLimit,    ///< simple upper limit for the value, e.g. p_T < 1 GeV/c
                     kAbsUpperLimit, ///< upper limit of the absolute value, e.g. |eta| < 0.8
                     kLowerLimit,    ///< simple lower limit for the value, e.g. p_T > 0.2 GeV/c
                     kAbsLowerLimit, ///< lower limit of the absolute value, e.g. |DCA_xyz| > 0.05 cm
                     kEqual          ///< values need to be equal, e.g. sign = 1
};

} // namespace femtoUniverseSelection

/// Simple class taking care of individual selections
/// \todo In principle all cuts that fulfill the getMinimalSelection are done implicitly and can be removed from the vector containing all cuts
/// \tparam selValDataType Data type used for the selection (float/int/bool/...)
/// \tparam selVariableDataType Data type used for the variable of the selection
template <class selValDataType, class selVariableDataType>
class FemtoUniverseSelection
{
 public:
  /// Default constructor
  FemtoUniverseSelection();

  /// Constructor
  /// \param selVal Value used for the selection
  /// \param selVar Variable used for the selection
  /// \param selType Type of selection to be employed
  FemtoUniverseSelection(selValDataType selVal, selVariableDataType selVar, femtoUniverseSelection::SelectionType selType)
    : mSelVal(selVal),
      mSelVar(selVar),
      mSelType(selType)
  {
  }

  /// Destructor
  virtual ~FemtoUniverseSelection() = default;

  /// Get the value used for the selection
  /// \return Value used for the selection
  selValDataType getSelectionValue() { return mSelVal; }

  /// Get the variable used for the selection
  /// \return variable used for the selection
  selVariableDataType getSelectionVariable() { return mSelVar; }

  /// Get the type of selection to be employed
  /// \return Type of selection to be employed
  femtoUniverseSelection::SelectionType getSelectionType() { return mSelType; }

  /// Check whether the selection is fulfilled or not
  /// \param observable Value of the variable to be checked
  /// \return Whether the selection is fulfilled or not
  bool isSelected(selValDataType observable)
  {
    switch (mSelType) {
      case (femtoUniverseSelection::SelectionType::kUpperLimit):
        return (observable < mSelVal);
      case (femtoUniverseSelection::SelectionType::kAbsUpperLimit):
        return (std::abs(observable) < mSelVal);
        break;
      case (femtoUniverseSelection::SelectionType::kLowerLimit):
        return (observable > mSelVal);
      case (femtoUniverseSelection::SelectionType::kAbsLowerLimit):
        return (std::abs(observable) > mSelVal);
        break;
      case (femtoUniverseSelection::SelectionType::kEqual):
        /// \todo can the comparison be done a bit nicer?
        return (std::abs(observable - mSelVal) < std::abs(mSelVal * 1e-6));
        break;
    }
    return false;
  }

  /// Check whether the selection is fulfilled or not and put together the bit-wise container for the systematic variations
  /// \tparam T Data type of the bit-wise container for the systematic variations
  /// \param observable Value of the variable to be checked
  /// \param cutContainer Bit-wise container for the systematic variations
  /// \param counter Position in the bit-wise container for the systematic variations to be modified
  template <typename T>
  void checkSelectionSetBit(selValDataType observable, T& cutContainer, size_t& counter)
  {
    /// If the selection is fulfilled the bit at the specified position (counter) within the bit-wise container is set to 1
    if (isSelected(observable)) {
      cutContainer |= 1UL << counter;
    }
    ++counter;
  }

 private:
  selValDataType mSelVal{0.f};                    ///< Value used for the selection
  selVariableDataType mSelVar;                    ///< Variable used for the selection
  femtoUniverseSelection::SelectionType mSelType; ///< Type of selection employed
};

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSESELECTION_H_
