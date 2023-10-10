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

/// \file FemtoDreamCutculator.h
/// \brief FemtoDreamCutculator - small class to match bit-wise encoding and
/// actual physics cuts \author Andi Mathis, TU München,
/// andreas.mathis@ph.tum.de \author Luca Barioglio, TU München,
/// luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_
#define PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_

#include "FemtoDreamSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "FemtoDreamV0Selection.h"
#include <bitset>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamCutculator
/// \brief Small class to match bit-wise encoding and actual physics cuts
class FemtoDreamCutculator
{
 public:
  /// Initializes boost ptree
  /// \param configFile Path to the dpl-config.json file from the
  /// femtodream-producer task
  void init(const char* configFile)
  {
    std::cout << "Welcome to the CutCulator!" << std::endl;

    boost::property_tree::ptree root;
    try {
      boost::property_tree::read_json(configFile, root);
    } catch (const boost::property_tree::ptree_error& e) {
      std::cout
        << "Failed to read JSON config file " << configFile << " ("
        << e.what() << ")" << std::endl;
    }

    // check the config file for all known producer task
    std::vector<const char*> ProducerTasks = {
      "femto-dream-producer-task", "femto-dream-producer-reduced-task"};
    for (auto& Producer : ProducerTasks) {
      if (root.count(Producer) > 0) {
        mConfigTree = root.get_child(Producer);
        std::cout << "Found " << Producer << " in " << configFile << std::endl;
        break;
      }
    }
  };

  /// Generic function that retrieves a given selection from the boost ptree and
  /// returns an std::vector in the proper format \param name Name of the
  /// selection in the dpl-config.json \return std::vector that can be directly
  /// passed to the FemtoDreamTrack/V0/../Selection
  std::vector<float> setSelection(std::string name)
  {
    try {
      boost::property_tree::ptree& selections = mConfigTree.get_child(name);
      boost::property_tree::ptree& selectionsValues =
        selections.get_child("values");
      std::vector<float> tmpVec;
      for (boost::property_tree::ptree::value_type& val : selectionsValues) {
        tmpVec.push_back(std::stof(val.second.data()));
      }
      return tmpVec;
    } catch (const boost::property_tree::ptree_error& e) {
      std::cout << "Selection " << name << " not available (" << e.what() << ")"
                << std::endl;
      return {};
    }
  }

  /// Specialization of the setSelection function for tracks

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setTrackSelection(femtoDreamTrackSelection::TrackSel obs,
                         femtoDreamSelection::SelectionType type,
                         const char* prefix)
  {
    auto tmpVec =
      setSelection(FemtoDreamTrackSelection::getSelectionName(obs, prefix));
    if (tmpVec.size() > 0) {
      mTrackSel.setSelection(tmpVec, obs, type);
    }
  }

  /// Automatically retrieves track selections from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setTrackSelectionFromFile(const char* prefix)
  {
    for (const auto& sel : mConfigTree) {
      std::string sel_name = sel.first;
      femtoDreamTrackSelection::TrackSel obs;
      if (sel_name.find(prefix) != std::string::npos) {
        int index = FemtoDreamTrackSelection::findSelectionIndex(
          std::string_view(sel_name), prefix);
        if (index >= 0) {
          obs = femtoDreamTrackSelection::TrackSel(index);
        } else {
          continue;
        }
        if (obs == femtoDreamTrackSelection::TrackSel::kPIDnSigmaMax)
          continue; // kPIDnSigmaMax is a special case
        setTrackSelection(obs, FemtoDreamTrackSelection::getSelectionType(obs),
                          prefix);
      }
    }
  }

  /// Automatically retrieves track PID from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setPIDSelectionFromFile(const char* prefix)
  {
    std::string PIDnodeName = std::string(prefix) + "PIDspecies";
    std::string PIDNsigmaNodeName = std::string(prefix) + "PIDnSigmaMax";
    try {
      boost::property_tree::ptree& pidNode = mConfigTree.get_child(PIDnodeName);
      boost::property_tree::ptree& pidValues = pidNode.get_child("values");
      for (auto& val : pidValues) {
        mPIDspecies.push_back(
          static_cast<o2::track::PID::ID>(std::stoi(val.second.data())));
      }
      boost::property_tree::ptree& pidNsigmaNode = mConfigTree.get_child(PIDNsigmaNodeName);
      boost::property_tree::ptree& pidNsigmaValues = pidNsigmaNode.get_child("values");
      for (auto& val : pidNsigmaValues) {
        mPIDValues.push_back(std::stof(val.second.data()));
      }
    } catch (const boost::property_tree::ptree_error& e) {
      std::cout << "PID selection not avalible for these skimmed data."
                << std::endl;
    }
  }

  /// Specialization of the setSelection function for V0

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setV0Selection(femtoDreamV0Selection::V0Sel obs,
                      femtoDreamSelection::SelectionType type,
                      const char* prefix)
  {
    auto tmpVec =
      setSelection(FemtoDreamV0Selection::getSelectionName(obs, prefix));
    if (tmpVec.size() > 0) {
      mV0Sel.setSelection(tmpVec, obs, type);
    }
  }

  /// Automatically retrieves V0 selections from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setV0SelectionFromFile(const char* prefix)
  {
    for (const auto& sel : mConfigTree) {
      std::string sel_name = sel.first;
      femtoDreamV0Selection::V0Sel obs;
      if (sel_name.find(prefix) != std::string::npos) {
        int index = FemtoDreamV0Selection::findSelectionIndex(
          std::string_view(sel_name), prefix);
        if (index >= 0) {
          obs = femtoDreamV0Selection::V0Sel(index);
        } else {
          continue;
        }
        setV0Selection(obs, FemtoDreamV0Selection::getSelectionType(obs),
                       prefix);
      }
    }
  }

  /// This function investigates a given selection criterion. The available
  /// options are displayed in the terminal and the bit-wise container is put
  /// together according to the user input \tparam T1 Selection class under
  /// investigation \param T2  Selection type under investigation \param output
  /// Bit-wise container for the systematic variations \param counter Current
  /// position in the bit-wise container to modify \tparam objectSelection
  /// Selection class under investigation (FemtoDreamTrack/V0/../Selection)
  /// \param selectionType Selection type under investigation, as defined in the
  /// selection class
  template <typename T1, typename T2>
  void checkForSelection(aod::femtodreamparticle::cutContainerType& output,
                         size_t& counter, T1 objectSelection, T2 selectionType,
                         bool SysChecks, float sign)
  {
    /// Output of the available selections and user input
    std::cout << "Selection: "
              << objectSelection.getSelectionHelper(selectionType) << " - (";
    auto selVec = objectSelection.getSelections(selectionType);
    for (auto selIt : selVec) {
      std::cout << selIt.getSelectionValue() << " ";
    }
    std::cout << ")" << std::endl
              << " > ";
    std::string in;
    std::vector<float> out;
    float input;

    if (SysChecks) {
      if (objectSelection.getSelectionHelper(selectionType) == std::string("Sign of the track")) {
        input = sign;
        std::cout << sign << std::endl;
      } else {
        // Seed the random number generator
        std::random_device rd;
        std::mt19937 rng(rd());
        // Select a random element
        std::uniform_int_distribution<int> uni(0, selVec.size() - 1);
        int randomIndex = uni(rng);
        input = selVec[randomIndex].getSelectionValue();
        std::cout << input << std::endl;
      }
    } else {
      if (selVec.size() == 1) {
        input = selVec[0].getSelectionValue();
        std::cout << input << std::endl;
      } else {
        std::cin >> in;
        input = std::stof(in);
      }
    }

    /// First we check whether the input is actually contained within the
    /// options
    bool inputSane = false;
    for (auto sel : selVec) {
      if (std::abs(sel.getSelectionValue() - input) <=
          std::abs(1.e-6 * input)) {
        inputSane = true;
      }
    }

    /// If the input is sane, the selection bit is put together
    if (inputSane) {
      int internal_index = 0;
      for (auto sel : selVec) {
        double signOffset;
        switch (sel.getSelectionType()) {
          case femtoDreamSelection::SelectionType::kEqual:
            signOffset = 0.;
            break;
          case (femtoDreamSelection::SelectionType::kLowerLimit):
          case (femtoDreamSelection::SelectionType::kAbsLowerLimit):
            signOffset = 1.;
            break;
          case (femtoDreamSelection::SelectionType::kUpperLimit):
          case (femtoDreamSelection::SelectionType::kAbsUpperLimit):
            signOffset = -1.;
            break;
        }

        /// for upper and lower limit we have to subtract/add an epsilon so that
        /// the cut is actually fulfilled
        if (sel.isSelected(input + signOffset * 1.e-6 * input)) {
          output |= 1UL << counter;
          for (int i = internal_index; i > 0; i--) {
            output &= ~(1UL << (counter - i));
          }
        }
        ++counter;
        ++internal_index;
      }
    } else {
      std::cout << "Choice " << in << " not recognized - repeating"
                << std::endl;
      checkForSelection(output, counter, objectSelection, selectionType, SysChecks, sign);
    }
  }

  /// This function iterates over all selection types of a given class and puts
  /// together the bit-wise container \tparam T1 Selection class under
  /// investigation \tparam objectSelection Selection class under investigation
  /// (FemtoDreamTrack/V0/../Selection) \return the full selection bit-wise
  /// container that will be put to the user task incorporating the user choice
  /// of selections
  template <typename T>
  aod::femtodreamparticle::cutContainerType iterateSelection(T objectSelection,
                                                             bool SysChecks, float sign)
  {
    aod::femtodreamparticle::cutContainerType output = 0;
    size_t counter = 0;
    auto selectionVariables = objectSelection.getSelectionVariables();
    for (auto selVarIt : selectionVariables) {
      checkForSelection(output, counter, objectSelection, selVarIt, SysChecks, sign);
    }
    return output;
  }

  /// This is the function called by the executable that then outputs the full
  /// selection bit-wise container incorporating the user choice of selections
  void analyseCuts(std::string choice, bool SysChecks = false, float sign = 1)
  {
    aod::femtodreamparticle::cutContainerType output = -1;
    if (choice == std::string("T")) {
      output = iterateSelection(mTrackSel, SysChecks, sign);
    } else if (choice == std::string("V")) {
      output = iterateSelection(mV0Sel, SysChecks, sign);
    } else {
      std::cout << "Option " << choice
                << " not recognized - available options are (T/V)" << std::endl;
      return;
    }
    std::bitset<8 * sizeof(aod::femtodreamparticle::cutContainerType)>
      bitOutput = output;
    std::cout << "+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "CutCulator has spoken - your selection bit is" << std::endl;
    std::cout << bitOutput << " (bitwise)" << std::endl;
    std::cout << output << " (number representation)" << std::endl;
    std::cout << "PID for these species is stored:" << std::endl;
    int index = 0;
    int randomIndex = 0;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, mPIDValues.size() - 1);
    for (auto id : mPIDspecies) {
      std::cout << o2::track::PID::getName(id) << " : " << index++ << std::endl;
      if (SysChecks) {
        // Seed the random number generator
        // Select a random element
        randomIndex = uni(rng);
        std::cout << "Nsigma TPC: " << mPIDValues[randomIndex] << std::endl;
        randomIndex = uni(rng);
        std::cout << "Nsigma TPCTOF: " << mPIDValues[randomIndex] << std::endl;
      }
    }
    std::cout << "+++++++++++++++++++++++++++++++++" << std::endl;
  }

 private:
  boost::property_tree::ptree
    mConfigTree; ///< the dpl-config.json buffered into a ptree
  FemtoDreamTrackSelection
    mTrackSel; ///< for setting up the bit-wise selection container for tracks
  FemtoDreamV0Selection
    mV0Sel; ///< for setting up the bit-wise selection container for V0s
  std::vector<o2::track::PID::ID>
    mPIDspecies; ///< list of particle species for which PID is stored
  std::vector<float>
    mPIDValues; ///< list of nsigma values for which PID is stored
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_ */
