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

/// \file FemtoUniverseCutculator.h
/// \brief FemtoUniverseCutculator - small class to match bit-wise encoding and actual physics cuts
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECUTCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECUTCULATOR_H_

#include <bitset>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"

namespace o2::analysis::femto_universe
{

/// \class FemtoUniverseCutculator
/// \brief Small class to match bit-wise encoding and actual physics cuts
class FemtoUniverseCutculator
{
 public:
  /// Initializes boost ptree
  /// \param configFile Path to the dpl-config.json file from the
  /// femtouniverse-producer task
  void init(const char* configFile)
  {
    LOGF(info, "Welcome to the CutCulator!");
    // std::cout << "Welcome to the CutCulator!" << std::endl;

    boost::property_tree::ptree root;
    try {
      boost::property_tree::read_json(configFile, root);
    } catch (const boost::property_tree::ptree_error& e) {
      // LOGF(fatal, "Failed to read JSON config file  %s (%s)", configFile, e.what());
      std::cout << "Failed to read JSON config file " << configFile << " (" << e.what() << ")" << std::endl;
    }

    // check the config file for all known producer task
    std::vector<const char*> producerTasks = {"femto-universe-producer-task"};
    for (const auto& Producer : producerTasks) {
      if (root.count(Producer) > 0) {
        mConfigTree = root.get_child(Producer);
        // LOGF(info, "Found %s in %s", Producer, configFile);
        std::cout << "Found " << Producer << " in " << configFile << std::endl;
        break;
      }
    }
  };

  /// Generic function that retrieves a given selection from the boost ptree and
  /// returns an std::vector in the proper format \param name Name of the
  /// selection in the dpl-config.json \return std::vector that can be directly
  /// passed to the FemtoUniverseTrack/V0/../Selection
  std::vector<float> setSelection(std::string name)
  {
    try {
      boost::property_tree::ptree& selections = mConfigTree.get_child(name);
      boost::property_tree::ptree& selectionsValues = selections.get_child("values");
      std::vector<float> tmpVec;
      for (boost::property_tree::ptree::value_type& val : selectionsValues) {
        tmpVec.push_back(std::stof(val.second.data()));
      }
      return tmpVec;
    } catch (const boost::property_tree::ptree_error& e) {
      // LOGF(fatal, "Selection %s not available (%s)", name, e.what());
      std::cout << "Selection " << name << " not available (" << e.what() << ")" << std::endl;
      return {};
    }
  }

  /// Specialization of the setSelection function for tracks

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setTrackSelection(femto_universe_track_selection::TrackSel obs,
                         femto_universe_selection::SelectionType type,
                         const char* prefix)
  {
    auto tmpVec = setSelection(FemtoUniverseTrackSelection::getSelectionName(obs, prefix));
    if (tmpVec.size() > 0) {
      mTrackSel.setSelection(tmpVec, obs, type);
    }
  }

  /// Automatically retrieves track selections from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setTrackSelectionFromFile(const char* prefix)
  {
    for (const auto& sel : mConfigTree) {
      std::string selName = sel.first;
      femto_universe_track_selection::TrackSel obs;
      if (selName.find(prefix) != std::string::npos) {
        int index = FemtoUniverseTrackSelection::findSelectionIndex(
          std::string_view(selName), prefix);
        if (index >= 0) {
          obs = femto_universe_track_selection::TrackSel(index);
        } else {
          continue;
        }
        if (obs == femto_universe_track_selection::TrackSel::kPIDnSigmaMax)
          continue; // kPIDnSigmaMax is a special case
        setTrackSelection(obs, FemtoUniverseTrackSelection::getSelectionType(obs),
                          prefix);
      }
    }
  }

  /// Automatically retrieves track PID from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setPIDSelectionFromFile(const char* prefix)
  {
    std::string mPIDnodeName = std::string(prefix) + "PIDspecies";
    std::string mPIDNsigmaNodeName = std::string(prefix) + "PIDnSigmaMax";
    try {
      boost::property_tree::ptree& pidNode = mConfigTree.get_child(mPIDnodeName);
      boost::property_tree::ptree& pidValues = pidNode.get_child("values");
      for (const auto& val : pidValues) {
        mPIDspecies.push_back(
          static_cast<o2::track::PID::ID>(std::stoi(val.second.data())));
      }
      boost::property_tree::ptree& pidNsigmaNode = mConfigTree.get_child(mPIDNsigmaNodeName);
      boost::property_tree::ptree& pidNsigmaValues = pidNsigmaNode.get_child("values");
      for (const auto& val : pidNsigmaValues) {
        mPIDValues.push_back(std::stof(val.second.data()));
      }
    } catch (const boost::property_tree::ptree_error& e) {
      // LOGF(fatal, "PID selection not avalible for these skimmed data.");
      std::cout << "PID selection not avalible for these skimmed data." << std::endl;
    }
  }

  /// Specialization of the setSelection function for V0

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setV0Selection(femto_universe_v0_selection::V0Sel obs,
                      femto_universe_selection::SelectionType type,
                      const char* prefix)
  {
    auto tmpVec =
      setSelection(FemtoUniverseV0Selection::getSelectionName(obs, prefix));
    if (tmpVec.size() > 0) {
      mV0Sel.setSelection(tmpVec, obs, type);
    }
  }

  /// Automatically retrieves V0 selections from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setV0SelectionFromFile(const char* prefix)
  {
    for (const auto& sel : mConfigTree) {
      std::string selName = sel.first;
      femto_universe_v0_selection::V0Sel obs;
      if (selName.find(prefix) != std::string::npos) {
        int index = FemtoUniverseV0Selection::findSelectionIndex(
          std::string_view(selName), prefix);
        if (index >= 0) {
          obs = femto_universe_v0_selection::V0Sel(index);
        } else {
          continue;
        }
        setV0Selection(obs, FemtoUniverseV0Selection::getSelectionType(obs),
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
  /// Selection class under investigation (FemtoUniverseTrack/V0/../Selection)
  /// \param selectionType Selection type under investigation, as defined in the
  /// selection class
  template <typename T1, typename T2>
  void checkForSelection(aod::femtouniverseparticle::CutContainerType& output,
                         size_t& counter, T1 objectSelection, T2 selectionType,
                         bool SysChecks, float sign)
  {
    /// Output of the available selections and user input
    std::cout << "Selection: " << objectSelection.getSelectionHelper(selectionType) << " - (";
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
      int internalIndex = 0;
      for (auto sel : selVec) {
        double signOffset;
        switch (sel.getSelectionType()) {
          case femto_universe_selection::SelectionType::kEqual:
            signOffset = 0.;
            break;
          case (femto_universe_selection::SelectionType::kLowerLimit):
          case (femto_universe_selection::SelectionType::kAbsLowerLimit):
            signOffset = 1.;
            break;
          case (femto_universe_selection::SelectionType::kUpperLimit):
          case (femto_universe_selection::SelectionType::kAbsUpperLimit):
            signOffset = -1.;
            break;
        }

        /// for upper and lower limit we have to subtract/add an epsilon so that
        /// the cut is actually fulfilled
        if (sel.isSelected(input + signOffset * 1.e-6 * input)) {
          output |= 1UL << counter;
          for (int i = internalIndex; i > 0; i--) {
            output &= ~(1UL << (counter - i));
          }
        }
        ++counter;
        ++internalIndex;
      }
    } else {
      // LOGF(info, "Choice %s not recognized - repeating", in);
      std::cout << "Choice " << in << " not recognized - repeating" << std::endl;
      checkForSelection(output, counter, objectSelection, selectionType, SysChecks, sign);
    }
  }

  /// This function iterates over all selection types of a given class and puts
  /// together the bit-wise container \tparam T1 Selection class under
  /// investigation \tparam objectSelection Selection class under investigation
  /// (FemtoUniverseTrack/V0/../Selection) \return the full selection bit-wise
  /// container that will be put to the user task incorporating the user choice
  /// of selections
  template <typename T>
  aod::femtouniverseparticle::CutContainerType iterateSelection(T objectSelection, bool SysChecks, float sign)
  {
    aod::femtouniverseparticle::CutContainerType output = 0;
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
    aod::femtouniverseparticle::CutContainerType output = -1;
    if (choice == std::string("T")) {
      output = iterateSelection(mTrackSel, SysChecks, sign);
    } else if (choice == std::string("V")) {
      output = iterateSelection(mV0Sel, SysChecks, sign);
    } else {
      // LOGF(info, "Option %s not recognized - available options are (T/V)", choice);
      std::cout << "Option " << choice << " not recognized - available options are (T/V)" << std::endl;
      return;
    }
    std::bitset<8 * sizeof(aod::femtouniverseparticle::CutContainerType)>
      bitOutput = output;
    // LOGF(info, "+++++++++++++++++++++++++++++++++");
    // LOGF(info, "CutCulator has spoken - your selection bit is");
    // LOGF(info, "%s (bitwise)", bitOutput);
    // LOGF(info, "%s (number representation)", output);
    // LOGF(info, "PID for these species is stored:");
    std::cout << "+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "CutCulator has spoken - your selection bit is" << std::endl;
    std::cout << bitOutput << " (bitwise)" << std::endl;
    std::cout << output << " (number representation)" << std::endl;
    std::cout << "PID for these species is stored:" << std::endl;
    int index = 0;
    for (auto id : mPIDspecies) {
      // LOGF(info, "%s : %d", o2::track::PID::getName(id), index++);
      std::cout << o2::track::PID::getName(id) << " : " << index++ << std::endl;
      if (SysChecks) {
        // Seed the random number generator
        std::random_device rd;
        std::mt19937 rng(rd());
        // Select a random element
        std::uniform_int_distribution<int> uni(0, mPIDValues.size() - 1);
        int randomIndex = uni(rng);
        std::cout << "Nsigma: " << mPIDValues[randomIndex] << std::endl;
      }
    }
  }

 private:
  boost::property_tree::ptree
    mConfigTree; ///< the dpl-config.json buffered into a ptree
  FemtoUniverseTrackSelection
    mTrackSel; ///< for setting up the bit-wise selection container for tracks
  FemtoUniverseV0Selection
    mV0Sel; ///< for setting up the bit-wise selection container for V0s
  std::vector<o2::track::PID::ID>
    mPIDspecies; ///< list of particle species for which PID is stored
  std::vector<float>
    mPIDValues; ///< list of nsigma values for which PID is stored
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECUTCULATOR_H_ */
