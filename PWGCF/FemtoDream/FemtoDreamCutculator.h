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

/// \file FemtoDreamCutculator.h
/// \brief FemtoDreamCutculator - small class to match bit-wise encoding and actual physics cuts
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_
#define PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_

#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "FemtoDreamSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "FemtoDreamV0Selection.h"

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamCutculator
/// \brief Small class to match bit-wise encoding and actual physics cuts
class FemtoDreamCutculator
{
 public:
  /// Initializes boost ptree
  /// \param configFile Path to the dpl-config.json file from the femtodream-producer task
  void init(const char* configFile)
  {
    LOG(info) << "Welcome to the CutCulator!";

    boost::property_tree::ptree root;
    try {
      boost::property_tree::read_json(configFile, root);
    } catch (const boost::property_tree::ptree_error& e) {
      LOG(fatal) << "Failed to read JSON config file " << configFile << " (" << e.what() << ")";
    }

    // check the config file for all known producer task
    std::vector<const char*> ProducerTasks = {"femto-dream-producer-task", "femto-dream-producer-reduced-task"};
    for (auto& Producer : ProducerTasks) {
      if (root.count(Producer) > 0) {
        mConfigTree = root.get_child(Producer);
        LOG(info) << "Found " << Producer << " in " << configFile;
        break;
      }
    }
  };

  /// Generic function that retrieves a given selection from the boost ptree and returns an std::vector in the proper format
  /// \param name Name of the selection in the dpl-config.json
  /// \return std::vector that can be directly passed to the FemtoDreamTrack/V0/../Selection
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
      LOG(warning) << "Selection " << name << " not available (" << e.what() << ")";
      return {};
    }
  }

  /// Specialization of the setSelection function for tracks

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setTrackSelection(femtoDreamTrackSelection::TrackSel obs, femtoDreamSelection::SelectionType type, const char* prefix)
  {
    auto tmpVec = setSelection(FemtoDreamTrackSelection::getSelectionName(obs, prefix));
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
        int index = FemtoDreamTrackSelection::findSelectionIndex(std::string_view(sel_name), prefix);
        if (index >= 0) {
          obs = femtoDreamTrackSelection::TrackSel(index);
        } else {
          continue;
        }
        if (obs == femtoDreamTrackSelection::TrackSel::kPIDnSigmaMax)
          continue; // kPIDnSigmaMax is a special case
        setTrackSelection(obs, FemtoDreamTrackSelection::getSelectionType(obs), prefix);
      }
    }
  }

  /// Automatically retrieves track PID from the dpl-config.json
  /// \param prefix Prefix which is added to the name of the Configurable
  void setPIDSelectionFromFile(const char* prefix)
  {
    std::string PIDnodeName = std::string(prefix) + "species";
    try {
      boost::property_tree::ptree& pidNode = mConfigTree.get_child(PIDnodeName);
      boost::property_tree::ptree& pidValues = pidNode.get_child("values");
      for (auto& val : pidValues) {
        mPIDspecies.push_back(static_cast<o2::track::PID::ID>(std::stoi(val.second.data())));
      }
    } catch (const boost::property_tree::ptree_error& e) {
      LOG(info) << "PID selection not avalible for these skimmed data.";
    }
  }

  /// Specialization of the setSelection function for V0

  /// The selection passed to the function is retrieved from the dpl-config.json
  /// \param obs Observable of the track selection
  /// \param type Type of the track selection
  /// \param prefix Prefix which is added to the name of the Configurable
  void setV0Selection(femtoDreamV0Selection::V0Sel obs, femtoDreamSelection::SelectionType type, const char* prefix)
  {
    auto tmpVec = setSelection(FemtoDreamV0Selection::getSelectionName(obs, prefix));
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
        int index = FemtoDreamV0Selection::findSelectionIndex(std::string_view(sel_name), prefix);
        if (index >= 0) {
          obs = femtoDreamV0Selection::V0Sel(index);
        } else {
          continue;
        }
        setV0Selection(obs, FemtoDreamV0Selection::getSelectionType(obs), prefix);
      }
    }
  }

  /// This function investigates a given selection criterion. The available options are displayed in the terminal and the bit-wise container is put together according to the user input
  /// \tparam T1 Selection class under investigation
  /// \param T2  Selection type under investigation
  /// \param output Bit-wise container for the systematic variations
  /// \param counter Current position in the bit-wise container to modify
  /// \tparam objectSelection Selection class under investigation (FemtoDreamTrack/V0/../Selection)
  /// \param selectionType Selection type under investigation, as defined in the selection class
  template <typename T1, typename T2>
  void checkForSelection(aod::femtodreamparticle::cutContainerType& output, size_t& counter, T1 objectSelection, T2 selectionType)
  {
    /// Output of the available selections and user input
    std::cout << "Selection: " << objectSelection.getSelectionHelper(selectionType) << " - (";
    auto selVec = objectSelection.getSelections(selectionType);
    for (auto selIt : selVec) {
      std::cout << selIt.getSelectionValue() << " ";
    }
    std::cout << ")\n > ";
    std::string in;
    std::cin >> in;
    const float input = std::stof(in);

    /// First we check whether the input is actually contained within the options
    bool inputSane = false;
    for (auto sel : selVec) {
      if (std::abs(sel.getSelectionValue() - input) <= std::abs(1.e-6 * input)) {
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

        /// for upper and lower limit we have to subtract/add an epsilon so that the cut is actually fulfilled
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
      std::cout << "Choice " << in << " not recognized - repeating\n";
      checkForSelection(output, counter, objectSelection, selectionType);
    }
  }

  /// This function iterates over all selection types of a given class and puts together the bit-wise container
  /// \tparam T1 Selection class under investigation
  /// \tparam objectSelection Selection class under investigation (FemtoDreamTrack/V0/../Selection)
  /// \return the full selection bit-wise container that will be put to the user task incorporating the user choice of selections
  template <typename T>
  aod::femtodreamparticle::cutContainerType iterateSelection(T objectSelection)
  {
    aod::femtodreamparticle::cutContainerType output = 0;
    size_t counter = 0;
    auto selectionVariables = objectSelection.getSelectionVariables();
    for (auto selVarIt : selectionVariables) {
      checkForSelection(output, counter, objectSelection, selVarIt);
    }
    return output;
  }

  /// This is the function called by the executable that then outputs the full selection bit-wise container incorporating the user choice of selections
  void analyseCuts()
  {
    std::cout << "Do you want to work with tracks/v0/cascade (T/V/C)?\n";
    std::cout << " > ";
    std::string in;
    std::cin >> in;
    aod::femtodreamparticle::cutContainerType output = -1;
    if (in.compare("T") == 0) {
      output = iterateSelection(mTrackSel);
    } else if (in.compare("V") == 0) {
      output = iterateSelection(mV0Sel);
    } else if (in.compare("C") == 0) {
      // output =  iterateSelection(mCascadeSel);
    } else {
      std::cout << "Option " << in << " not recognized - available options are (T/V/C) \n";
      analyseCuts();
    }
    std::bitset<8 * sizeof(aod::femtodreamparticle::cutContainerType)> bitOutput = output;
    std::cout << "+++++++++++++++++++++++++++++++++\n";
    std::cout << "CutCulator has spoken - your selection bit is\n";
    std::cout << bitOutput << " (bitwise)\n";
    std::cout << output << " (number representation)\n";
    std::cout << "PID for these species is stored:\n";
    int index = 0;
    for (auto id : mPIDspecies) {
      std::cout << o2::track::PID::getName(id) << " : " << index++ << std::endl;
    }
  }

 private:
  boost::property_tree::ptree mConfigTree;     ///< the dpl-config.json buffered into a ptree
  FemtoDreamTrackSelection mTrackSel;          ///< for setting up the bit-wise selection container for tracks
  FemtoDreamV0Selection mV0Sel;                ///< for setting up the bit-wise selection container for V0s
  std::vector<o2::track::PID::ID> mPIDspecies; ///< list of particle species for which PID is stored
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMCUTCULATOR_H_ */
