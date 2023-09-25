#!/usr/bin/env python3

"""
Generates a C++ header with PDG codes and particle masses.
Author: Vít Kučera <vit.kucera@cern.ch>
"""

import ROOT
from ctypes import c_bool
from enum import Enum

class PdgROOT(Enum):
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

def mass(code):
    if abs(code) == Pdg.kXiCCPlusPlus.value:
        return 3.62155  # https://pdg.lbl.gov/ (2021)
    if abs(code) == Pdg.kOmegaC0.value:
        return 2.6952  # https://pdg.lbl.gov/ (2022)
    success = c_bool(True)
    return ROOT.o2.O2DatabasePDG.Mass(code, success)

type = "double"

str_enum = "enum Code {\n"
str_mass = ""
for c in Pdg:
    str_enum += f"  {c.name} = {c.value},\n"
    str_mass += f"constexpr {type} Mass{c.name[1:]} = {mass(c.value)};\n"
str_enum += "};\n"
print(str_enum)
print(str_mass)

str_mass_root = ""
for c in PdgROOT:
    str_mass_root += f"constexpr {type} Mass{c.name[1:]} = {mass(c.value)};\n"
print(str_mass_root)
