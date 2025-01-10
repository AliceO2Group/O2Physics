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

/// \file femtoUniverseCutCulator.cxx
/// \brief Executable that encodes physical selection criteria in a bit-wise selection
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include <filesystem>
#include <iostream>
#include <random>
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCutculator.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

using namespace o2::analysis::femto_universe;

/// The function takes the path to the dpl-config.json as a argument and the
/// does a Q&A session for the user to find the appropriate selection criteria
/// for the analysis task
int main(int /*argc*/, char* argv[])
{
  std::string configFileName(argv[1]);
  std::filesystem::path configFile{configFileName};

  if (std::filesystem::exists(configFile)) {
    FemtoUniverseCutculator cut;
    cut.init(argv[1]);

    std::cout
      << "Do you want to work with tracks or V0s (T/V)? >";
    std::string choice;
    std::cin >> choice;

    if (choice == std::string("T")) {
      cut.setTrackSelectionFromFile("ConfTrk");
      cut.setPIDSelectionFromFile("ConfTrk");
    } else if (choice == std::string("V")) {
      std::cout << "Do you want to select V0s or one of its children (V/T)? >";
      std::cin >> choice;
      cut.setV0SelectionFromFile("ConfV0");
      cut.setTrackSelectionFromFile("ConfChild");
      cut.setPIDSelectionFromFile("ConfChild");
    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }
    /// \todo factor out the pid here
    /// cut.setTrackSelection(femto_universe_track_selection::kPIDnSigmaMax,
    /// femto_universe_selection::kAbsUpperLimit, "ConfTrk");

    std::cout << "Do you want to manually select cuts or create systematic "
                 "variations(M/V)? >";
    std::string manual;
    std::cin >> manual;

    if (manual == std::string("M")) {
      cut.analyseCuts(choice);
    } else if (manual == std::string("V")) {
      std::ofstream out("CutCulator.txt");
      std::streambuf* coutbuf = std::cout.rdbuf(); // save old buf
      std::cout.rdbuf(out.rdbuf());                // redirect std::cout to out.txt!
      for (int i = 0; i < 20; i++) {
        cut.analyseCuts(choice, true, 1);
      }
      std::cout.rdbuf(coutbuf); // reset to standard output again
    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }

  } else {
    std::cout << "The configuration file " << configFileName
              << " could not be found.";
  }

  return 0;
}
