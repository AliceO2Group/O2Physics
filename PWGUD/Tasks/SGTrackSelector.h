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
//
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#ifndef PWGUD_TASKS_SGTRACKSELECTOR_H_
#define PWGUD_TASKS_SGTRACKSELECTOR_H_

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "PWGUD/Core/RLhelper.h"
#include <TString.h>
#include "TLorentzVector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396
#define mkaon 0.4937
#define mproton 0.9383
template <typename T>
int trackselector(const T& track)
{
  TLorentzVector a;
  a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
  if (std::abs(track.dcaZ()) > 2.)
    return 0;
  if (std::abs(track.dcaXY()) > .0105 + .035 / pow(a.Pt(), 1.1))
    return 0;
  if (track.tpcChi2NCl() > 4)
    return 0;
  if (track.tpcNClsFindable() < 70)
    return 0;
  if (track.itsChi2NCl() > 36)
    return 0;
  return 1;
}
template <typename T>
int trackpid(const T& track, bool use_tof)
{
  int pid = 0;
  float pi, ka, pr;
  float tpi, tka, tpr;
  pi = std::abs(track.tpcNSigmaPi());
  ka = std::abs(track.tpcNSigmaKa());
  pr = std::abs(track.tpcNSigmaPr());
  if (pi < 1. && pi < ka && pi < pr)
    pid = 1;
  else if (ka < 1. && ka < pi && ka < pr)
    pid = 2;
  else if (pr < 1. && pr < pi && pr < ka)
    pid = 3;
  if (use_tof && track.tofChi2() > -1) {
    tpi = std::abs(track.tofNSigmaPi());
    tka = std::abs(track.tofNSigmaKa());
    tpr = std::abs(track.tofNSigmaPr());
    if (std::sqrt(pi * pi + tpi * tpi) < 2 && std::sqrt(pi * pi + tpi * tpi) < std::sqrt(ka * ka + tka * tka) && std::sqrt(pi * pi + tpi * tpi) < std::sqrt(pr * pr + tpr * tpr))
      pid = 1;
    else if (std::sqrt(ka * ka + tka * tka) < 2 && std::sqrt(pi * pi + tpi * tpi) > std::sqrt(ka * ka + tka * tka) && std::sqrt(ka * ka + tka * tka) < std::sqrt(pr * pr + tpr * tpr))
      pid = 2;
    else if (std::sqrt(pr * pr + tpr * tpr) < 2 && std::sqrt(pr * pr + tpr * tpr) < std::sqrt(ka * ka + tka * tka) && std::sqrt(pi * pi + tpi * tpi) > std::sqrt(pr * pr + tpr * tpr))
      pid = 3;
  }
  return pid;
}

#endif // PWGUD_TASKS_SGTRACKSELECTOR_H_
