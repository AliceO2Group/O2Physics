// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamCutCulator.cxx
/// \brief Executable that encodes physical selection criteria in a bit-wise selection
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include <filesystem>
#include <iostream>
#include "FemtoDreamCutculator.h"
#include "FemtoDreamSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2::analysis::femtoDream;

/// The function takes the path to the dpl-config.json as a argument and the
/// does a Q&A session for the user to find the appropriate selection criteria
/// for the analysis task
int main(int argc, char* argv[])
{
  std::string configFileName(argv[1]);
  std::filesystem::path configFile{configFileName};

  if (std::filesystem::exists(configFile)) {
    FemtoDreamCutculator cut;
    cut.init(argv[1]);

    LOG(info) << "Do you want to work with tracks or V0s (T/V)?";
    std::string choice;
    std::cin >> choice;

    if (choice == std::string("T")) {
      cut.setTrackSelectionFromFile("ConfTrk");
      cut.setPIDSelectionFromFile("ConfPIDTrk");
    } else if (choice == std::string("V")) {
      LOG(info) << "Do you want to select V0s or one of its children (V/T)?";
      std::cin >> choice;
      cut.setV0SelectionFromFile("ConfV0");
      cut.setTrackSelectionFromFile("ConfChild");
      cut.setPIDSelectionFromFile("ConfPIDChild");
    } else {
      LOG(info) << "Option not recognized. Break...";
      return 1;
    }

    /// \todo factor out the pid here
    // cut.setTrackSelection(femtoDreamTrackSelection::kPIDnSigmaMax,
    // femtoDreamSelection::kAbsUpperLimit, "ConfTrk");
    cut.analyseCuts(choice);

  } else {
    LOG(info) << "The configuration file " << configFileName
              << " could not be found.";
  }
}
