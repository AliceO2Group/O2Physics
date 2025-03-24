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

/// \file FemtoDreamSelection.h
/// \brief FemtoDreamSelection - small generic class to do selections
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMSELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMSELECTION_H_

#include <cmath>
#include "Framework/HistogramRegistry.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femtoDream
{

namespace femtoDreamSelection
{
/// Type of selection to be employed
enum SelectionType { kUpperLimit,    ///< simple upper limit for the value, e.g. p_T < 1 GeV/c
                     kAbsUpperLimit, ///< upper limit of the absolute value, e.g. |eta| < 0.8
                     kLowerLimit,    ///< simple lower limit for the value, e.g. p_T > 0.2 GeV/c
                     kAbsLowerLimit, ///< lower limit of the absolute value, e.g. |DCA_xyz| > 0.05 cm
                     kEqual          ///< values need to be equal, e.g. sign = 1
};

} // namespace femtoDreamSelection

/// Simple class taking care of individual selections
/// \todo In principle all cuts that fulfill the getMinimalSelection are done implicitly and can be removed from the vector containing all cuts
/// \tparam selValDataType Data type used for the selection (float/int/bool/...)
/// \tparam selVariableDataType Data type used for the variable of the selection
template <class selValDataType, class selVariableDataType>
class FemtoDreamSelection
{
 public:
  /// Default constructor
  FemtoDreamSelection();

  /// Constructor
  /// \param selVal Value used for the selection
  /// \param selVar Variable used for the selection
  /// \param selType Type of selection to be employed
  FemtoDreamSelection(selValDataType selVal, selVariableDataType selVar, femtoDreamSelection::SelectionType selType)
    : mSelVal(selVal),
      mSelVar(selVar),
      mSelType(selType)
  {
  }

  /// Destructor
  virtual ~FemtoDreamSelection() = default;

  /// Get the value used for the selection
  /// \return Value used for the selection
  selValDataType getSelectionValue() { return mSelVal; }

  /// Get the variable used for the selection
  /// \return variable used for the selection
  selVariableDataType getSelectionVariable() { return mSelVar; }

  /// Get the type of selection to be employed
  /// \return Type of selection to be employed
  femtoDreamSelection::SelectionType getSelectionType() { return mSelType; }

  /// Check whether the selection is fulfilled or not
  /// \param observable Value of the variable to be checked
  /// \return Whether the selection is fulfilled or not
  bool isSelected(selValDataType observable)
  {
    switch (mSelType) {
      case (femtoDreamSelection::SelectionType::kUpperLimit):
        return (observable <= mSelVal);
      case (femtoDreamSelection::SelectionType::kAbsUpperLimit):
        return (std::fabs(observable) <= mSelVal);
        break;
      case (femtoDreamSelection::SelectionType::kLowerLimit):
        return (observable >= mSelVal);
      case (femtoDreamSelection::SelectionType::kAbsLowerLimit):
        return (std::fabs(observable) >= mSelVal);
        break;
      case (femtoDreamSelection::SelectionType::kEqual):
        /// \todo can the comparison be done a bit nicer?
        return (std::fabs(observable - mSelVal) < std::abs(mSelVal * 1e-6));
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
  void checkSelectionSetBit(selValDataType observable, T& cutContainer, size_t& counter, HistogramRegistry* registry)
  {
    /// If the selection is fulfilled the bit at the specified position (counter) within the bit-wise container is set to 1
    if (isSelected(observable)) {
      cutContainer |= 1UL << counter;
      if (registry) {
        registry->fill(HIST("AnalysisQA/CutCounter"), 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType));
      }
    } else {
      if (registry) {
        registry->fill(HIST("AnalysisQA/CutCounter"), counter);
      }
    }
    ++counter;
  }

  template <typename T>
  void checkSelectionSetBitPID(selValDataType observable, T& cutContainer)
  {
    /// If the selection is fulfilled the bit at the specified position (counter) within the bit-wise container is set to 1
    if (isSelected(observable)) {
      cutContainer |= 1UL;
    } else {
      cutContainer |= 0UL;
    }
    cutContainer <<= 1;
  }

 private:
  selValDataType mSelVal{0.f};                 ///< Value used for the selection
  selVariableDataType mSelVar;                 ///< Variable used for the selection
  femtoDreamSelection::SelectionType mSelType; ///< Type of selection employed
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMSELECTION_H_
