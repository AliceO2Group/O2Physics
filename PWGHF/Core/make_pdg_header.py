#!/usr/bin/env python3

"""
Generates the body of a C++ header with PDG codes and particle masses.
Author: Vít Kučera <vit.kucera@cern.ch>
"""

import ROOT  # pylint: disable=import-error
from ctypes import c_bool
from enum import Enum

# Enum of PDG_t particles
class PdgROOT(Enum):
    kDown = ROOT.kDown
    kDownBar = ROOT.kDownBar
    kUp = ROOT.kUp
    kUpBar = ROOT.kUpBar
    kStrange = ROOT.kStrange
    kStrangeBar = ROOT.kStrangeBar
    kCharm = ROOT.kCharm
    kCharmBar = ROOT.kCharmBar
    kBottom = ROOT.kBottom
    kBottomBar = ROOT.kBottomBar
    kTop = ROOT.kTop
    kTopBar = ROOT.kTopBar
    kGluon = ROOT.kGluon
    kElectron = ROOT.kElectron
    kPositron = ROOT.kPositron
    kNuE = ROOT.kNuE
    kNuEBar = ROOT.kNuEBar
    kMuonMinus = ROOT.kMuonMinus
    kMuonPlus = ROOT.kMuonPlus
    kNuMu = ROOT.kNuMu
    kNuMuBar = ROOT.kNuMuBar
    kTauMinus = ROOT.kTauMinus
    kTauPlus = ROOT.kTauPlus
    kNuTau = ROOT.kNuTau
    kNuTauBar = ROOT.kNuTauBar
    kGamma = ROOT.kGamma
    kZ0 = ROOT.kZ0
    kWPlus = ROOT.kWPlus
    kWMinus = ROOT.kWMinus
    kPi0 = ROOT.kPi0
    kK0Long = ROOT.kK0Long
    kPiPlus = ROOT.kPiPlus
    kPiMinus = ROOT.kPiMinus
    kProton = ROOT.kProton
    kProtonBar = ROOT.kProtonBar
    kNeutron = ROOT.kNeutron
    kNeutronBar = ROOT.kNeutronBar
    kK0Short = ROOT.kK0Short
    kK0 = ROOT.kK0
    kK0Bar = ROOT.kK0Bar
    kKPlus = ROOT.kKPlus
    kKMinus = ROOT.kKMinus
    kLambda0 = ROOT.kLambda0
    kLambda0Bar = ROOT.kLambda0Bar
    kLambda1520 = ROOT.kLambda1520
    kSigmaMinus = ROOT.kSigmaMinus
    kSigmaBarPlus = ROOT.kSigmaBarPlus
    kSigmaPlus = ROOT.kSigmaPlus
    kSigmaBarMinus = ROOT.kSigmaBarMinus
    kSigma0 = ROOT.kSigma0
    kSigma0Bar = ROOT.kSigma0Bar
    kXiMinus = ROOT.kXiMinus
    kXiPlusBar = ROOT.kXiPlusBar
    kOmegaMinus = ROOT.kOmegaMinus
    kOmegaPlusBar = ROOT.kOmegaPlusBar

# Enum of additional particles
class Pdg(Enum):
    kB0 = 511
    kB0Bar = -511
    kBPlus = 521
    kBS = 531
    kBSBar = -531
    kD0 = 421
    kD0Bar = -421
    kDMinus = -411
    kDPlus = 411
    kDS = 431
    kDSBar = -431
    kDStar = 413
    kChiC1 = 20443
    kJPsi = 443
    kLambdaB0 = 5122
    kLambdaCPlus = 4122
    kOmegaC0 = 4332
    kPhi = 333
    kSigmaC0 = 4112
    kSigmaCPlusPlus = 4222
    kX3872 = 9920443
    kXiCCPlusPlus = 4422
    kXiCPlus = 4232
    kXiCZero = 4132

"""
Returns particle mass from o2::O2DatabasePDG except for special cases.
"""
def mass(code):
    # Special cases (present in TDatabasePDG but with wrong values)
    # Missing particles should be added in O2DatabasePDG.h.
    if abs(code) == Pdg.kXiCCPlusPlus.value:
        return 3.62155  # PDG 2021
    if abs(code) == Pdg.kOmegaC0.value:
        return 2.69520  # PDG 2022
    # Default case
    success = c_bool(True)
    return ROOT.o2.O2DatabasePDG.Mass(code, success)

"""
Returns a C++ declaration of a particle mass constant.
"""
def declare_mass(pdg, type="double") -> str:
    return f"constexpr {type} Mass{pdg.name[1:]} = {mass(pdg.value)};\n"

# Start of enum declarations of additional particles
str_enum_head = """/// \\brief Declarations of named PDG codes of particles missing in ROOT PDG_t
/// \\note Follow kCamelCase naming convention
/// \\link https://root.cern/doc/master/TPDGCode_8h.html
enum Code {
"""
# End of enum declarations of additional particles
str_enum_foot = "};\n"
# Documentation string for mass declarations of additional particles
str_mass_o2_head = """/// \\brief Declarations of masses for additional particles
"""
# Documentation string for mass declarations of PDG_t particles
str_mass_root_head = """/// \\brief Declarations of masses for particles in ROOT PDG_t
"""

# Additional particles
str_enum = str_enum_head
str_mass_o2 = str_mass_o2_head
for c in Pdg:
    str_enum += f"  {c.name} = {c.value},\n"
    str_mass_o2 += declare_mass(c)
str_enum = str_enum[:-2] + "\n"  # Remove the last comma.
str_enum += str_enum_foot

# PDG_t particles
str_mass_root = str_mass_root_head
for d in PdgROOT:
    str_mass_root += declare_mass(d)

# Header body
str_header = "\n".join([str_enum, str_mass_o2, str_mass_root])
print(str_header)
