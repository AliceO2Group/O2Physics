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

/// \file FemtoWorldObjectSelection.h
/// \brief FemtoWorldObjectSelection - Parent class of all selections
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef FEMTOWORLDOBJECTSELECTION_H_
#define FEMTOWORLDOBJECTSELECTION_H_

#include "PWGCF/FemtoWorld/Core/FemtoWorldSelection.h"

#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2::framework;

namespace o2::analysis
{
namespace femtoWorld
{

/// \class FemtoWorldObjectSelection
/// \brief Cut class to contain and execute all cuts applied to tracks
/// \todo In principle all cuts that fulfill the getMinimalSelection are done implicitly and can be removed from the vector containing all cuts
/// \tparam selValDataType Data type used for the selection (float/int/bool/...)
/// \tparam selVariable Variable used for the selection
template <class selValDataType, class selVariable>
class FemtoWorldObjectSelection
{
 public:
  /// Destructor
  virtual ~FemtoWorldObjectSelection() = default;

  /// The selection criteria employed in the child class are written to a histogram
  /// \tparam part Type of the particle, used as a prefix for the folder in the QAResults.root
  template <o2::aod::femtoworldparticle::ParticleType part>
  void fillSelectionHistogram()
  {
    int nBins = mSelections.size();
    mHistogramRegistry->add((static_cast<std::string>(o2::aod::femtoworldparticle::ParticleTypeName[part]) + "/cuthist").c_str(), "; Cut; Value", kTH1F, {{nBins, 0, static_cast<double>(nBins)}});
    auto hist = mHistogramRegistry->get<TH1>(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/cuthist"));
    for (size_t i = 0; i < mSelections.size(); ++i) {
      hist->GetXaxis()->SetBinLabel(i + 1, Form("%u", mSelections.at(i).getSelectionVariable()));
      hist->SetBinContent(i + 1, mSelections.at(i).getSelectionValue());
    }
  }

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \tparam T Type of the configurable passed to the function
  /// \param selVals o2 configurable containing the values employed for the selection
  /// \param selVar Variable to be employed for the selection
  /// \param selType Type of the selection to be employed
  template <typename T>
  void setSelection(T& selVals, selVariable selVar, femtoWorldSelection::SelectionType selType)
  {
    std::vector<selValDataType> tmpSelVals = selVals; // necessary due to some features of the Configurable
    std::vector<FemtoWorldSelection<selValDataType, selVariable>> tempVec;
    for (const selValDataType selVal : tmpSelVals) {
      tempVec.push_back(FemtoWorldSelection<selValDataType, selVariable>(selVal, selVar, selType));
    }
    setSelection(tempVec);
  }

  /// Pass an std::vector of selection values to the selection class
  /// \param sels std::vector containing FemtoWorldSelections
  void setSelection(std::vector<FemtoWorldSelection<selValDataType, selVariable>>& sels)
  {
    /// First the selection is sorted so that the most open cuts are conducted first
    switch (sels.at(0).getSelectionType()) {
      case (femtoWorldSelection::SelectionType::kUpperLimit):
      case (femtoWorldSelection::SelectionType::kAbsUpperLimit):
        std::sort(sels.begin(), sels.end(), [](FemtoWorldSelection<selValDataType, selVariable> a, FemtoWorldSelection<selValDataType, selVariable> b) {
          return a.getSelectionValue() > b.getSelectionValue();
        });
        break;
      case (femtoWorldSelection::SelectionType::kLowerLimit):
      case (femtoWorldSelection::SelectionType::kAbsLowerLimit):
      case (femtoWorldSelection::SelectionType::kEqual):
        std::sort(sels.begin(), sels.end(), [](FemtoWorldSelection<selValDataType, selVariable> a, FemtoWorldSelection<selValDataType, selVariable> b) {
          return a.getSelectionValue() < b.getSelectionValue();
        });
        break;
    }

    /// Then, the sorted selections are added to the overall container of cuts
    for (const auto& sel : sels) {
      mSelections.push_back(sel);
    }
  }

  /// Retrieve the most open selection of a given selection variable
  /// \param selVar Selection variable under consideration
  /// \param selType Type of the selection variable
  /// \return The most open selection of the selection variable given to the class
  selValDataType getMinimalSelection(selVariable selVar, femtoWorldSelection::SelectionType selType)
  {
    selValDataType minimalSel{};
    switch (selType) {
      case (femtoWorldSelection::SelectionType::kUpperLimit):
      case (femtoWorldSelection::SelectionType::kAbsUpperLimit):
        minimalSel = -999.e9;
        break;
      case (femtoWorldSelection::SelectionType::kLowerLimit):
      case (femtoWorldSelection::SelectionType::kAbsLowerLimit):
      case (femtoWorldSelection::SelectionType::kEqual):
        minimalSel = 999.e9;
        break;
    }

    for (auto sel : mSelections) {
      if (sel.getSelectionVariable() == selVar) {
        switch (sel.getSelectionType()) {
          case (femtoWorldSelection::SelectionType::kUpperLimit):
          case (femtoWorldSelection::SelectionType::kAbsUpperLimit):
            if (minimalSel < sel.getSelectionValue()) {
              minimalSel = sel.getSelectionValue();
            }
            break;
          case (femtoWorldSelection::SelectionType::kLowerLimit):
          case (femtoWorldSelection::SelectionType::kAbsLowerLimit):
          case (femtoWorldSelection::SelectionType::kEqual):
            if (minimalSel > sel.getSelectionValue()) {
              minimalSel = sel.getSelectionValue();
            }
            break;
        }
      }
    }
    return minimalSel;
  }

  /// The total number of different selections
  /// \return Total number of selections
  size_t getNSelections()
  {
    return mSelections.size();
  }

  /// The number of selection of an individual variable
  /// \param selVar Selection variable under consideration
  /// \return Number of selection of the individual variable
  size_t getNSelections(selVariable selVar)
  {
    return getSelections(selVar).size();
  }

  /// Obtain the selections of an individual variable
  /// \param selVar Selection variable under consideration
  /// \return All selections of the individual variable
  std::vector<FemtoWorldSelection<selValDataType, selVariable>> getSelections(selVariable selVar)
  {
    std::vector<FemtoWorldSelection<selValDataType, selVariable>> selValVec;
    for (auto it : mSelections) {
      if (it.getSelectionVariable() == selVar) {
        selValVec.push_back(it);
      }
    }
    return selValVec;
  }

  /// Retrieve all the different selection variables
  /// \return std::vector containing all the different selection variables
  std::vector<selVariable> getSelectionVariables()
  {
    std::vector<selVariable> selVarVec;
    for (auto it : mSelections) {
      auto selVar = it.getSelectionVariable();
      if (std::none_of(selVarVec.begin(), selVarVec.end(), [selVar](selVariable a) { return a == selVar; })) {
        selVarVec.push_back(selVar);
      }
    }
    return selVarVec;
  }

 protected:
  HistogramRegistry* mHistogramRegistry;                                     ///< For QA output
  std::vector<FemtoWorldSelection<selValDataType, selVariable>> mSelections; ///< Vector containing all selections
};

} // namespace femtoWorld
} // namespace o2::analysis

#endif /* FEMTOWORLDOBJECTSELECTION_H_ */
