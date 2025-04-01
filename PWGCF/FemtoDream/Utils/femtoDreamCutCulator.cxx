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

/// \file femtoDreamCutCulator.cxx
/// \brief Executable that encodes physical selection criteria in a bit-wise
/// selection \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include <iostream>
#include <random>
#include <string>
#include "PWGCF/FemtoDream/Utils/femtoDreamCutCulator.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2::analysis::femtoDream;

/// The function takes the path to the dpl-config.json as a argument and the
/// does a Q&A session for the user to find the appropriate selection criteria
/// for the analysis task
int main(int /*argc*/, char* argv[])
{
  std::string configFileName(argv[1]);
  std::ifstream configFile(configFileName);

  if (configFile.is_open()) {
    FemtoDreamCutculator cut;
    cut.init(argv[1]);

    std::cout
      << "Do you want to work with tracks or V0s or Cascades (T/V/C)? >";
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
    } else if (choice == std::string("C")) {
      std::cout << "Do you want to select cascades, V0-Daughter tracks of the cascades or the Bachelor track (C/V/B)? >";
      std::cin >> choice;
      if (choice == std::string("C")) {
        cut.setCascadeSelectionFromFile("ConfCascade");
        choice = "C";
      } else if (choice == std::string("V")) {
        cut.setTrackSelectionFromFile("ConfCascV0Child");
        cut.setPIDSelectionFromFile("ConfCascV0Child");
        choice = "T";
      } else if (choice == std::string("B")) {
        cut.setTrackSelectionFromFile("ConfCascBachelor");
        cut.setPIDSelectionFromFile("ConfCascBachelor");
        choice = "T";
      } else {
        std::cout << "Option not recognized. Break...";
        return 2;
      }
    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }
    /// \todo factor out the pid here
    /// cut.setTrackSelection(femtoDreamTrackSelection::kPIDnSigmaMax,
    /// femtoDreamSelection::kAbsUpperLimit, "ConfTrk");

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
      for (int i = 0; i < 1000; i++) {
        cut.analyseCuts(choice, true, 1);
      }
      std::cout.rdbuf(coutbuf); // reset to standard output again
    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }

  } else {
    std::cout << "The configuration file " << configFileName
              << " could not be found or could not be opened.";
    return 1;
  }

  return 0;
}
