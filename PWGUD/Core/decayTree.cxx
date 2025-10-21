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

#include "decayTree.h"

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace rapidjson;

// -----------------------------------------------------------------------------
// pidSelector holds an array of pidcut
// pidcut
// double[8]
//  0:  pid to apply
//  1:  detector: 1 - TPC, 2 - TOF
//  2:  cut type:
//     1: pt and nSigma within limits
//    -1: nSigma out of limits within pt range
//     2: pt and detector signal  within limits
//    -2: detector signal out of limits within pt range
//  3:  How to apply cut:
//     0: not active
//     1: if information available
//     2: return false if information not available
//  4:  pT min
//  5:  pT max
//  6:  signal/nSigma min
//  7:  signal/nSigma max
// pidSelector
pidSelector::pidSelector(std::vector<std::vector<double>>& pidcuts)
{
  fpidCuts = pidcuts;
}

void pidSelector::clear()
{
  fpidCuts.clear();
}

int pidSelector::pid2ind(int pid)
{
  switch (std::abs(pid)) {
    case 11: // electron
      return 0;
    case 211: // pion
      return 1;
    case 13: // muon
      return 2;
    case 321: // kaon
      return 3;
    case 2212: // proton
      return 4;
    default: // unknown
      return -1.;
  }
};

void pidSelector::Print()
{
  LOGF(info, "      PID cuts");
  for (const auto& cut : fpidCuts) {
    LOGF(info, "        [ %.0f %.0f %.0f %.0f %.0f %.2f %.2f %.2f ]", cut[0], cut[1], cut[2], cut[3], cut[4], cut[5], cut[6], cut[7]);
  }
}

// -----------------------------------------------------------------------------
// angleCut
angleCut::angleCut(std::pair<std::string, std::string> rnames, double angleMin, double angleMax)
{
  fRnames = rnames;
  fAngleMin = angleMin;
  fAngleMax = angleMax;
}

void angleCut::Print()
{
  LOGF(info, "        %s ^ %s: %f : %f", fRnames.first, fRnames.second, fAngleMin, fAngleMax);
}

// -----------------------------------------------------------------------------
// Resonance
resonance::resonance()
{
  init();
}

// reset to default values
void resonance::init()
{
  // initialisations
  fisFinal = false;
  fCounter = -1;
  fName = "";
  fStatus = 0;
  fPID = 0;
  fPIDfun = 0;
  fdetectorHits = std::vector<int>{-1, -1, -1, -1};
  fParents.clear();
  fDaughters.clear();

  // initialize mass, pT, eta range
  setMassRange(0., 100.);
  setPtRange(0., 100.);
  setEtaRange(-10., 10.);
  setNcltpcRange(0, 200);
  setChi2ncltpcRange(0., 100.);
  setDCAxyzMax(100., 100.);

  // histogram axes
  fnmassBins = 350;
  fmassHistMin = 0.0;
  fmassHistMax = 3.5;
  fnmomBins = 300;
  fmomHistMin = 0.0;
  fmomHistMax = 3.0;

  // invariant mass
  fIVM = TLorentzVector(0., 0., 0., 0.);
  fCharge = 0;

  // pid cuts
  fpidSelector.clear();
  fangleCuts.clear();
}

void resonance::reset()
{
  fStatus = 0;
}

// check mass, pt, eta range and charge
void resonance::updateStatus()
{
  // IVM has to be computed
  if (fStatus == 0) {
    return;
  }

  // check mass, pt, and eta range
  fStatus = 2;
  if (fIVM.M() < fmassMin || fIVM.M() > fmassMax) {
    return;
  }
  if (fIVM.Perp() < fptMin || fIVM.Perp() > fptMax) {
    return;
  }
  if (fIVM.Eta() < fetaMin || fIVM.Eta() > fetaMax) {
    return;
  }

  // good candidate
  fStatus = 3;
}

void resonance::Print()
{
  if (fisFinal) {
    // final
    LOGF(info, "    %s : %d", fName, fPID);
    LOGF(info, "      status: %d", fStatus);
    LOGF(info, "      final: %d", fCounter);
    LOGF(info, "      nCluster TPC: %d : %d", fncltpcMin, fncltpcMax);
    LOGF(info, "      chi2 per cluster TPC: %f : %f", fchi2ncltpcMin, fchi2ncltpcMax);
    LOGF(info, "      maximum dca_XY : %f", fdcaxyMax);
    LOGF(info, "      maximum dca_Z: %f", fdcazMax);
    LOGF(info, "      parents");
    for (const auto& parent : fParents) {
      LOGF(info, "        %s", parent);
    }
  } else {
    // resonance
    LOGF(info, "    %s", fName);
    LOGF(info, "      status: %d", fStatus);
    LOGF(info, "      parents");
    for (const auto& parent : fParents) {
      LOGF(info, "        %s", parent);
    }
    LOGF(info, "      daughters");
    for (const auto& daugh : fDaughters) {
      LOGF(info, "        %s", daugh);
    }
  }
  LOGF(info, "      mass range: %f : %f", fmassMin, fmassMax);
  LOGF(info, "      pt range: %f : %f", fptMin, fptMax);
  LOGF(info, "      eta range: %f : %f", fetaMin, fetaMax);
  if (fisFinal) {
    fpidSelector.Print();
  } else {
    LOGF(info, "      Angle cuts");
    for (const auto& anglecut : fangleCuts) {
      anglecut->Print();
    }
  }

  if (fStatus > 0) {
    LOGF(info, "        charge: %d", fCharge);
    LOGF(info, "        mass: %f", fIVM.M());
    LOGF(info, "        E: %f", fIVM.E());
    LOGF(info, "        pT: %f", fIVM.Perp());
    LOGF(info, "        eta: %f", fIVM.Eta());
  }
  LOGF(info, "");
}

// -----------------------------------------------------------------------------
// decayTree
decayTree::decayTree()
{
  fPDG = TDatabasePDG::Instance();
}

bool decayTree::init(std::string const& parFile, o2::framework::HistogramRegistry& registry)
{
  // initialisation of constants
  fccs = {"ULS", "LS"};
  fdets = {"TPC", "TOF"};
  fparts = {"el", "pi", "mu", "ka", "pr"};

  // open the decayTree file
  FILE* fjson = fopen(parFile.c_str(), "r");
  if (!fjson) {
    LOGF(error, "Could not open parameter file %s", parFile);
    return false;
  }

  // create streamer
  char readBuffer[65536];
  FileReadStream jsonStream(fjson, readBuffer, sizeof(readBuffer));

  // parse the json file
  Document jsonDocument;
  jsonDocument.ParseStream(jsonStream);

  // is it a proper json document?
  if (jsonDocument.HasParseError()) {
    LOGF(error, "Check the parameter file! There is a problem with the format!");
    return false;
  }

  // check for decayTree
  const char* itemName = "decayTree";
  if (!jsonDocument.HasMember(itemName)) {
    LOGF(error, "Check the parameter file! Item %s is missing!", itemName);
    return false;
  }
  const Value& decTree = jsonDocument[itemName];

  // finals and resonance parameters
  std::string name;
  int pid;
  int pidfun;
  std::vector<std::string> daughters;
  std::vector<int> charges;
  std::vector<std::vector<double>> vcuts;
  std::vector<angleCut*> vanglecuts;

  itemName = "finals";
  if (!decTree.HasMember(itemName)) {
    LOGF(info, "No %s defined!", itemName);
    return false;
  }
  if (!decTree[itemName].IsArray()) {
    LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
    return false;
  }

  // loop over finals
  fnFinals = 0;
  auto fins = decTree[itemName].GetArray();
  for (const auto& fin : fins) {
    if (!fin.IsObject()) {
      LOGF(error, "Check the parameter file! %s must be objects!", itemName);
      return false;
    }

    // create new finals (resonance)
    resonance* newRes = new resonance();
    newRes->setisFinal();

    // check for name
    itemName = "name";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsString()) {
        name = fin[itemName].GetString();
        newRes->setName(name);
      } else {
        LOGF(error, "Check the parameter file! %s must be a string!", itemName);
        delete newRes;
        return false;
      }
    } else {
      LOGF(error, "Check the parameter file! Finals must have a %s!", itemName);
      delete newRes;
      return false;
    }

    // check for pid
    itemName = "pid";
    if (fin.HasMember(itemName)) {
      pid = fin[itemName].GetInt();
      newRes->setPID(pid);
    } else {
      LOGF(error, "Check the parameter file! Finals must have a %s!", itemName);
      delete newRes;
      return false;
    }

    // check for ptrange
    itemName = "ptrange";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto lims = fin[itemName].GetArray();
        if (lims.Size() == 2) {
          newRes->setPtRange(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for etarange
    itemName = "etarange";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto lims = fin[itemName].GetArray();
        if (lims.Size() == 2) {
          newRes->setEtaRange(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for detectorhits
    // 0: ITS
    // 1: TPC
    // 2: TRD
    // 3: TOF
    itemName = "detectorhits";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto hits = fin[itemName].GetArray();
        if (hits.Size() == 4) {
          newRes->setDetectorHits(hits[0].GetInt(), hits[1].GetInt(), hits[2].GetInt(), hits[3].GetInt());
        } else {
          LOGF(error, "Check the parameter file! %s must have four elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for ncltpcrange
    itemName = "ncltpcrange";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto lims = fin[itemName].GetArray();
        if (lims.Size() == 2) {
          newRes->setNcltpcRange(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for chi2ncltpcrange
    itemName = "chi2ncltpcrange";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto lims = fin[itemName].GetArray();
        if (lims.Size() == 2) {
          newRes->setChi2ncltpcRange(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for dcaxyzmax
    itemName = "dcaxyzmax";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto lims = fin[itemName].GetArray();
        if (lims.Size() == 2) {
          newRes->setDCAxyzMax(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // pid cuts
    vcuts.clear();
    itemName = "pidcuts";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto pidcuts = fin[itemName].GetArray();
        for (const auto& pidcut : pidcuts) {
          if (pidcut.HasMember("pidcut")) {
            if (pidcut["pidcut"].IsArray()) {
              auto vals = pidcut["pidcut"].GetArray();
              if (vals.Size() == 8) {
                std::vector<double> vs;
                for (const auto& val : vals) {
                  vs.push_back(val.GetFloat());
                }
                vcuts.push_back(vs);
              } else {
                LOGF(error, "Check the parameter file! A pidcut has %d members", vals.Size());
                delete newRes;
                return false;
              }
            } else {
              LOGF(error, "Check the parameter file! Item pidcuts['pidcut'] must be an array!");
              delete newRes;
              return false;
            }
          }
        }
        auto pidsel = pidSelector(vcuts);
        newRes->setPIDSelector(pidsel);
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for massaxis
    itemName = "massaxis";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto max = fin[itemName].GetArray();
        if (max.Size() == 3) {
          newRes->setMassHistAxis(max[0].GetInt(), max[1].GetFloat(), max[2].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have three elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // check for momaxis
    itemName = "momaxis";
    if (fin.HasMember(itemName)) {
      if (fin[itemName].IsArray()) {
        auto momax = fin[itemName].GetArray();
        if (momax.Size() == 3) {
          newRes->setMomHistAxis(momax[0].GetInt(), momax[1].GetFloat(), momax[2].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! %s must have three elements!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
        delete newRes;
        return false;
      }
    }

    // update fResonances
    fnFinals++;
    newRes->setCounter(fnFinals - 1);
    fResonances.push_back(newRes);
  }

  itemName = "ulsstates";
  if (decTree.HasMember(itemName)) {
    LOGF(info, "No %s defined!", itemName);
    if (!decTree[itemName].IsArray()) {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }
    auto ulsstates = decTree[itemName].GetArray();
    for (const auto& obj : ulsstates) {
      if (obj.IsArray()) {
        auto ulsstate = obj.GetArray();
        if (static_cast<int>(ulsstate.Size()) != fnFinals) {
          LOGF(error, "%s with wrong number of elements!", itemName);
          return false;
        }
        charges.clear();
        for (const auto& ch : ulsstate) {
          charges.push_back(ch.GetInt());
        }
        LOGF(info, "adding charge state %d", chargeState(charges));
        fULSstates.push_back(chargeState(charges));
      } else {
        LOGF(error, "Check the parameter file! Elements of item %s must be arrays!", itemName);
        return false;
      }
    }
  }

  itemName = "lsstates";
  if (decTree.HasMember(itemName)) {
    LOGF(info, "No %s defined!", itemName);
    if (!decTree[itemName].IsArray()) {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }
    auto lsstates = decTree[itemName].GetArray();
    for (const auto& obj : lsstates) {
      if (obj.IsArray()) {
        auto lsstate = obj.GetArray();
        if (static_cast<int>(lsstate.Size()) != fnFinals) {
          LOGF(error, "%s with wrong number of elements!", itemName);
          return false;
        }
        charges.clear();
        for (const auto& ch : lsstate) {
          charges.push_back(ch.GetInt());
        }
        fLSstates.push_back(chargeState(charges));
      } else {
        LOGF(error, "Check the parameter file! Elements of item %s must be arrays!", itemName);
        return false;
      }
    }
  }

  // update permutations
  permutations(fnFinals, fPermutations);

  itemName = "resonances";
  if (!decTree.HasMember(itemName)) {
    LOGF(info, "No resonances defined!");
  } else {
    if (!decTree[itemName].IsArray()) {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }

    // loop over resonances
    auto ress = decTree[itemName].GetArray();
    for (const auto& res : ress) {
      if (!res.IsObject()) {
        LOGF(error, "Check the parameter file! %s must be objects!", itemName);
        return false;
      }

      // create new resonance
      resonance* newRes = new resonance();

      // check for name
      itemName = "name";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsString()) {
          name = res[itemName].GetString();
          newRes->setName(name);
        } else {
          LOGF(error, "Check the parameter file! %s must be a string!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Resonances must have a name!");
        delete newRes;
        return false;
      }

      // check for daughters
      itemName = "daughters";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto daughs = res[itemName].GetArray();
          if (daughs.Size() >= 2) {
            daughters.clear();
            for (const auto& daugh : daughs) {
              daughters.push_back(daugh.GetString());
            }
            newRes->setDaughters(daughters);
          } else {
            LOGF(error, "Check the parameter file! %s must have at least two elements!", itemName);
            delete newRes;
            return false;
          }
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      } else {
        LOGF(error, "Check the parameter file! Resonances must have daughters!");
        delete newRes;
        return false;
      }

      // check for massrange
      itemName = "massrange";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto lims = res[itemName].GetArray();
          if (lims.Size() != 2) {
            LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
            delete newRes;
            return false;
          }
          newRes->setMassRange(lims[0].GetFloat(), lims[1].GetFloat());
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // check for ptrange
      itemName = "ptrange";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto lims = res[itemName].GetArray();
          if (lims.Size() == 2) {
            newRes->setPtRange(lims[0].GetFloat(), lims[1].GetFloat());
          } else {
            LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
            delete newRes;
            return false;
          }
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // check for etarange
      itemName = "etarange";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto lims = res[itemName].GetArray();
          if (lims.Size() == 2) {
            newRes->setEtaRange(lims[0].GetFloat(), lims[1].GetFloat());
          } else {
            LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
            delete newRes;
            return false;
          }
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // check for pidfun
      itemName = "pidfun";
      if (res.HasMember(itemName)) {
        pidfun = res[itemName].GetInt();
        newRes->setPIDFun(pidfun);
      }

      // angle cuts
      vanglecuts.clear();
      itemName = "anglecuts";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto anglecuts = res[itemName].GetArray();
          for (const auto& anglecut : anglecuts) {
            if (anglecut.HasMember("pair")) {
              if (anglecut["pair"].IsArray()) {
                auto daughters = anglecut["pair"].GetArray();
                if (daughters.Size() == 2) {
                  std::pair<std::string, std::string> pair{daughters[0].GetString(), daughters[1].GetString()};
                  if (anglecut.HasMember("anglerange")) {
                    if (anglecut["anglerange"].IsArray()) {
                      auto arange = anglecut["anglerange"].GetArray();
                      if (arange.Size() == 2) {
                        vanglecuts.push_back(new angleCut(pair, arange[0].GetFloat(), arange[1].GetFloat()));
                      } else {
                        LOGF(error, "Check the parameter file! An anglecut['anglerange'] has %d members", arange.Size());
                        delete newRes;
                        return false;
                      }
                    } else {
                      LOGF(error, "Check the parameter file! Item anglerange must be an array!");
                      delete newRes;
                      return false;
                    }
                  } else {
                    LOGF(error, "Check the parameter file! An anglecut must have an anglerange!");
                    delete newRes;
                    return false;
                  }
                } else {
                  LOGF(error, "Check the parameter file! A anglecut['pair'] has %d members", daughters.Size());
                  delete newRes;
                  return false;
                }
              } else {
                LOGF(error, "Check the parameter file! Item anglecut must be an array!");
                delete newRes;
                return false;
              }
            } else {
              LOGF(error, "Check the parameter file! An anglecut must have a pair!");
              delete newRes;
              return false;
            }
          }
          newRes->setAngleCuts(vanglecuts);
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // check for massaxis
      itemName = "massaxis";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto max = res[itemName].GetArray();
          if (max.Size() == 3) {
            newRes->setMassHistAxis(max[0].GetInt(), max[1].GetFloat(), max[2].GetFloat());
          } else {
            LOGF(error, "Check the parameter file! %s must have three elements!", itemName);
            delete newRes;
            return false;
          }
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // check for momaxis
      itemName = "momaxis";
      if (res.HasMember(itemName)) {
        if (res[itemName].IsArray()) {
          auto momax = res[itemName].GetArray();
          if (momax.Size() == 3) {
            newRes->setMomHistAxis(momax[0].GetInt(), momax[1].GetFloat(), momax[2].GetFloat());
          } else {
            LOGF(error, "Check the parameter file! %s must have three elements!", itemName);
            delete newRes;
            return false;
          }
        } else {
          LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
          delete newRes;
          return false;
        }
      }

      // use the items to create a new resonance
      fResonances.push_back(newRes);
    }
  }

  // update the parents information
  updateParents();

  // check for eventCuts
  // set default values
  fnTracksMin = fnFinals;
  fnTracksMax = fnFinals;
  frgtwtofMin = 0.0;
  fdBCMin = 0;
  fdBCMax = 0;
  fFITvetos = {0, 1, 1, 0, 0};

  itemName = "eventCuts";
  if (!jsonDocument.HasMember(itemName)) {
    return true;
  }
  const Value& evset = jsonDocument[itemName];

  // check for ntrackrange
  itemName = "ntrackrange";
  if (evset.HasMember(itemName)) {
    if (evset[itemName].IsArray()) {
      auto lims = evset[itemName].GetArray();
      if (lims.Size() == 2) {
        fnTracksMin = lims[0].GetFloat();
        fnTracksMax = lims[1].GetFloat();
      } else {
        LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
        return false;
      }
    } else {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }
  }

  // check for rgtwtofmin
  itemName = "rgtrwtofmin";
  if (evset.HasMember(itemName)) {
    frgtwtofMin = evset[itemName].GetFloat();
  }

  // check for dBCrange
  itemName = "dBCrange";
  if (evset.HasMember(itemName)) {
    if (evset[itemName].IsArray()) {
      auto lims = evset[itemName].GetArray();
      if (lims.Size() == 2) {
        fdBCMin = lims[0].GetInt();
        fdBCMax = lims[1].GetInt();
      } else {
        LOGF(error, "Check the parameter file! %s must have two elements!", itemName);
        return false;
      }
    } else {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }
  }

  // check for FITvetos
  itemName = "FITvetos";
  if (evset.HasMember(itemName)) {
    if (evset[itemName].IsArray()) {
      auto vetoes = evset[itemName].GetArray();
      if (vetoes.Size() == 5) {
        fFITvetos.clear();
        for (auto i = 0; i < 5; i++) {
          fFITvetos.push_back(vetoes[i].GetInt());
        }
      } else {
        LOGF(error, "Check the parameter file! %s must have five elements!", itemName);
        return false;
      }
    } else {
      LOGF(error, "Check the parameter file! Item %s must be an array!", itemName);
      return false;
    }
  }

  // clean up
  fclose(fjson);

  // create the histograms
  createHistograms(registry);

  return true;
}

void decayTree::updateParents()
{
  // reset parents
  for (const auto& res : fResonances) {
    res->clearParents();
  }

  for (const auto& res : fResonances) {
    for (const auto& daughName : res->getDaughters()) {
      auto daugh = getResonance(daughName);
      daugh->addParent(res->name());
    }
  }
}

void decayTree::reset()
{
  fStatus = 0;
  fChargeState = -1;
  for (const auto& res : fResonances) {
    res->reset();
  }
}

// decayTree status
//  0: unset
//  1: not accepted
//  2: ULS accepted
//  3: LS accepted
void decayTree::updateStatus()
{
  bool isULS = true;
  bool isLS = true;
  fStatus = 1;
  for (const auto& res : fResonances) {
    if (res->status() < 3) {
      return;
    }

    // check the charge state
    updateChargeState();
    if (std::find(fULSstates.begin(), fULSstates.end(), fChargeState) == fULSstates.end()) {
      isULS = false;
    }
    if (std::find(fLSstates.begin(), fLSstates.end(), fChargeState) == fLSstates.end()) {
      isLS = false;
    }
  }
  if (isULS) {
    fStatus = 2;
  } else if (isLS) {
    fStatus = 3;
  }
}

void decayTree::Print()
{
  LOGF(info, "eventCuts");
  LOGF(info, "  ntrkRange: %d : %d", fnTracksMin, fnTracksMax);
  LOGF(info, "  dBCRange: %d : %d", fdBCMin, fdBCMax);
  LOGF(info, "  FITVetoes: [ %d %d %d %d %d ]", fFITvetos[0], fFITvetos[1], fFITvetos[2], fFITvetos[3], fFITvetos[4]);
  LOGF(info, "  ULS states");
  for (const auto& chstat : fULSstates) {
    LOGF(info, "    %d", chstat);
  }
  LOGF(info, "  LS states");
  for (const auto& chstat : fLSstates) {
    LOGF(info, "    %d", chstat);
  }
  LOGF(info, "");
  LOGF(info, "decayTree");
  LOGF(info, "  nResonances: %d", fResonances.size());
  LOGF(info, "  nFinals: %d", fnFinals);
  LOGF(info, "  Resonances");
  for (const auto& res : fResonances) {
    res->Print();
  }
}

resonance* decayTree::getResonance(std::string name)
{
  for (const auto& res : fResonances) {
    if (res->name() == name) {
      return res;
    }
  }
  LOGF(error, "A resonance %s does not exists!", name);
  return new resonance(); // is needed to satisfy return type
}

resonance* decayTree::getFinal(int counter)
{
  for (const auto& res : fResonances) {
    if (res->counter() == counter) {
      return res;
    }
  }
  LOGF(error, "The final %d does not exists!", counter);
  return new resonance(); // is needed to satisfy return type
}

std::vector<resonance*> decayTree::getFinals(resonance* res)
{
  std::vector<resonance*> resFinals;

  for (const auto& d1Name : res->getDaughters()) {
    auto d1 = getResonance(d1Name);
    if (d1->isFinal()) {
      resFinals.push_back(d1);
    } else {
      for (const auto& d2Name : d1->getDaughters()) {
        auto d2 = getResonance(d2Name);
        for (const auto& d3 : getFinals(d2)) {
          resFinals.push_back(d3);
        }
      }
    }
  }
  return resFinals;
}

// apply anglecuts
void decayTree::checkAngles()
{
  // loop over resonances
  for (const auto& res : fResonances) {
    auto anglecuts = res->getAngleCuts();
    // loop over angle cuts
    for (const auto& anglecut : anglecuts) {
      auto rnames = anglecut->rNames();
      auto anglerange = anglecut->angleRange();

      // compute angle between two resonances
      auto lv1 = getResonance(rnames.first)->IVM();
      auto lv2 = getResonance(rnames.second)->IVM();
      auto ang = lv1.Angle(lv2.Vect());

      // apply cut
      if (ang < anglerange.first || ang > anglerange.second) {
        res->setStatus(2);
        break;
      }
    }
  }
}

// compute charge state
int decayTree::chargeState(std::vector<int> chs)
{
  int chargeState = -1;
  if (static_cast<int>(chs.size()) == fnFinals) {
    // loop over elements of chargestate
    chargeState = 0;
    for (auto ind = 0; ind < fnFinals; ind++) {
      chargeState += (chs[ind] > 0) * std::pow(2, ind);
    }
  }
  return chargeState;
}
void decayTree::updateChargeState()
{
  fChargeState = 0;

  // loop over all finals
  for (auto ind = 0; ind < fnFinals; ind++) {
    fChargeState += (getFinal(ind)->charge() > 0) * std::pow(2, ind);
  }
}

std::size_t decayTree::combHash(std::vector<int>& comb)
{
  // comb contains indices to tracks
  // the combination hash is created from a string which is comprised of a sorted list of track indices

  // sort comb
  std::map<int, int> m_unsorted;
  for (size_t cnt = 0; cnt < comb.size(); cnt++) {
    m_unsorted.insert(std::pair<int, int>(comb[cnt], cnt));
  }
  std::vector<int> v_sorted(comb.size());
  partial_sort_copy(begin(comb), end(comb), begin(v_sorted), end(v_sorted));

  // create the hash
  std::hash<std::string> hasher;
  std::string hashstr{""};
  for (size_t cnt = 0; cnt < comb.size(); cnt++) {
    hashstr += std::to_string(v_sorted[cnt]);
  }

  return hasher(hashstr);
}

// find all permutations of n0 elements
void decayTree::permutations(std::vector<int>& ref, int n0, int np, std::vector<std::vector<int>>& perms)
{
  // create local reference
  auto ref2u = ref;

  // loop over np-1 rotations of last np elements of ref
  for (auto ii = 0; ii < np; ii++) {

    // create a new permutation
    // copy first n0-np elements from ref
    // then rotate last np elements of ref
    std::vector<int> perm(n0, 0);
    for (auto ii = 0; ii < n0 - np; ii++) {
      perm[ii] = ref2u[ii];
    }
    for (auto ii = n0 - np + 1; ii < n0; ii++) {
      perm[ii - 1] = ref2u[ii];
    }
    perm[n0 - 1] = ref2u[n0 - np];

    // add new permutation to the list of permuutations
    if (ii < (np - 1)) {
      perms.push_back(perm);
    }

    // if np>2 then do permutation of next level
    // use the new combination as reference
    if (np > 2) {
      auto newnp = np - 1;
      permutations(perm, n0, newnp, perms);
    }

    // update reference
    ref2u = perm;
  }
}

// find all permutations of n0 elements
int decayTree::permutations(int n0, std::vector<std::vector<int>>& perms)
{
  // initialize with first trivial combination
  perms.clear();
  if (n0 == 0) {
    return 0;
  }

  std::vector<int> ref(n0, 0);
  for (auto ii = 0; ii < n0; ii++) {
    ref[ii] = ii;
  }
  perms.push_back(ref);

  // iterate recursively
  permutations(ref, n0, n0, perms);

  return perms.size();
}

void decayTree::combinations(int n0, std::vector<int>& pool, int np, std::vector<int>& inds, int n,
                             std::vector<std::vector<int>>& combs)
{
  // loop over pool
  for (auto ii = 0; ii < n0 - n; ii++) {

    inds[n] = pool[ii];

    // if all inds are defined then print them out
    // else get next inds
    if (np == 1) {

      std::vector<int> comb(n + 1, 0);
      for (uint ii = 0; ii < inds.size(); ii++) {
        comb[ii] = inds[ii];
      }
      combs.push_back(comb);

    } else {

      auto n0new = n0 - ii;
      std::vector<int> newpool(n0new, 0);
      for (auto kk = 0; kk < n0new; kk++) {
        newpool[kk] = pool[kk + ii + 1];
      }

      auto npnew = np - 1;
      auto nnew = n + 1;
      combinations(n0new, newpool, npnew, inds, nnew, combs);
    }
  }
}

// find all possible selections of np out of n0
int decayTree::combinations(int n0, int np, std::vector<std::vector<int>>& combs)
{
  // initialisations
  combs.clear();
  if (n0 < np) {
    return 0;
  }

  std::vector<int> pool(n0, 0);
  for (auto ii = 0; ii < n0; ii++) {
    pool[ii] = ii;
  }
  std::vector<int> inds(np, 0);

  // iterate recursively
  combinations(n0, pool, np, inds, 0, combs);

  return combs.size();
}

std::vector<std::vector<int>> decayTree::combinations(int nPool)
{
  // all selections of fnFinals items from nPool elements
  std::vector<std::vector<int>> combs;
  combinations(nPool, fnFinals, combs);

  // permute the combinations
  std::vector<std::vector<int>> copes;
  for (const auto& comb : combs) {
    for (const auto& perm : fPermutations) {
      std::vector<int> cope(fnFinals, 0);
      for (auto jj = 0; jj < fnFinals; jj++) {
        cope[perm[jj]] = comb[jj];
      }
      copes.push_back(cope);
    }
  }

  return copes;
}

// -----------------------------------------------------------------------------
