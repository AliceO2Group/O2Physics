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

/// \file listMasses.h
/// \brief Generate list of masses from the PDG database.

#include <TPDGCode.h>

#include "SimulationDataFormat/O2DatabasePDG.h"
#include "CommonConstants/PhysicsConstants.h"

// #include "PWGHF/Core/PDG.h"
#include "PDG.h"

void listMasses()
{
bool success = false;
double mass = 0;

using Pdg = o2::O2DatabasePDG;
using namespace o2::analysis::pdg;
using namespace o2::constants::physics;

std::map<std::string, int> mapROOT = {
  {"Rootino", kRootino},
  {"Down", kDown},
  {"DownBar", kDownBar},
  {"Up", kUp},
  {"UpBar", kUpBar},
  {"Strange", kStrange},
  {"StrangeBar", kStrangeBar},
  {"Charm", kCharm},
  {"CharmBar", kCharmBar},
  {"Bottom", kBottom},
  {"BottomBar", kBottomBar},
  {"Top", kTop},
  {"TopBar", kTopBar},
  {"Gluon", kGluon},
  {"Pythia92", kPythia92},
  {"Dd1", kDd1},
  {"Dd1Bar", kDd1Bar},
  {"Ud0", kUd0},
  {"Ud0Bar", kUd0Bar},
  {"Ud1", kUd1},
  {"Ud1Bar", kUd1Bar},
  {"Uu1", kUu1},
  {"Uu1Bar", kUu1Bar},
  {"Sd0", kSd0},
  {"Sd0Bar", kSd0Bar},
  {"Sd1", kSd1},
  {"Sd1Bar", kSd1Bar},
  {"Su0", kSu0},
  {"Su0Bar", kSu0Bar},
  {"Su1", kSu1},
  {"Su1Bar", kSu1Bar},
  {"Searches0", kSearches0},
  {"Electron", kElectron},
  {"Positron", kPositron},
  {"NuE", kNuE},
  {"NuEBar", kNuEBar},
  {"MuonMinus", kMuonMinus},
  {"MuonPlus", kMuonPlus},
  {"NuMu", kNuMu},
  {"NuMuBar", kNuMuBar},
  {"TauMinus", kTauMinus},
  {"TauPlus", kTauPlus},
  {"NuTau", kNuTau},
  {"NuTauBar", kNuTauBar},
  {"Gamma", kGamma},
  {"Z0", kZ0},
  {"WPlus", kWPlus},
  {"WMinus", kWMinus},
  {"Pi0", kPi0},
  {"Rho770_0", kRho770_0},
  {"A2_1320_0", kA2_1320_0},
  {"Rho3_1690_0", kRho3_1690_0},
  {"K0Long", kK0Long},
  {"PiPlus", kPiPlus},
  {"PiMinus", kPiMinus},
  {"Rho770Plus", kRho770Plus},
  {"Rho770Minus", kRho770Minus},
  {"A2_1320Plus", kA2_1320Plus},
  {"Proton", kProton},
  {"ProtonBar", kProtonBar},
  {"Neutron", kNeutron},
  {"NeutronBar", kNeutronBar},
  {"K0Short", kK0Short},
  {"K0", kK0},
  {"K0Bar", kK0Bar},
  {"KPlus", kKPlus},
  {"KMinus", kKMinus},
  {"Lambda0", kLambda0},
  {"Lambda1520", kLambda1520},
  {"Lambda0Bar", kLambda0Bar},
  {"SigmaMinus", kSigmaMinus},
  {"SigmaBarPlus", kSigmaBarPlus},
  {"SigmaPlus", kSigmaPlus},
  {"SigmaBarMinus", kSigmaBarMinus},
  {"Sigma0", kSigma0},
  {"Sigma0Bar", kSigma0Bar},
  {"XiMinus", kXiMinus},
  {"XiPlusBar", kXiPlusBar},
  {"OmegaMinus", kOmegaMinus},
  {"OmegaPlusBar", kOmegaPlusBar}
};

std::map<std::string, int> map = {
  {"B0", kB0},
  {"B0Bar", kB0Bar},
  {"BPlus", kBPlus},
  {"BS", kBS},
  {"BSBar", kBSBar},
  {"D0", kD0},
  {"D0Bar", kD0Bar},
  {"DMinus", kDMinus},
  {"DPlus", kDPlus},
  {"DS", kDS},
  {"DSBar", kDSBar},
  {"DStar", kDStar},
  {"ChiC1", kChiC1},
  {"JPsi", kJPsi},
  {"LambdaB0", kLambdaB0},
  {"LambdaCPlus", kLambdaCPlus},
  {"OmegaC0", kOmegaC0},
  {"Phi", kPhi},
  {"SigmaC0", kSigmaC0},
  {"SigmaCPlusPlus", kSigmaCPlusPlus},
  {"X3872", kX3872},
  {"XiCCPlusPlus", kXiCCPlusPlus},
  {"XiCPlus", kXiCPlus},
  {"XiCZero", kXiCZero}
};



// for (const auto& [name, pdg] : map) {
for (const auto& [name, pdg] : mapROOT) {
  float mass = Pdg::Mass(pdg, success);
  Printf("constexpr float Mass%s = %g;", name.c_str(), mass);
}
  Printf("\nconstexpr float MassAlpha = %g;", 3.7273794 - MassAlpha);
}
