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

#ifndef PWGUD_CORE_SGTRACKSELECTOR_H_
#define PWGUD_CORE_SGTRACKSELECTOR_H_

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TVector3.h"
#include <TString.h>

#include <iostream>
#include <vector>

template <typename T>
int trackselector(const T& track, const std::vector<float>& params)
{
  // Ensure the params vector contains all the necessary parameters

  if (params.size() < 8) {
    throw std::runtime_error("Insufficient parameters provided");
  }
  TVector3 a;
  a.SetXYZ(track.px(), track.py(), track.pz());
  if (params[0] == 1 && !track.isPVContributor())
    return 0;
  if (std::abs(track.dcaZ()) > params[1])
    return 0;
  if (!params[2]) {
    if (std::abs(track.dcaXY()) > .0105 + .035 / pow(a.Pt(), 1.1))
      return 0;
  } else {
    if (std::abs(track.dcaXY()) > params[2])
      return 0;
  }
  if (track.tpcChi2NCl() > params[3])
    return 0;
  if (track.tpcNClsFindable() < params[4])
    return 0;
  if (track.itsChi2NCl() > params[5])
    return 0;
  if (std::abs(a.Eta()) > params[6])
    return 0;
  if (a.Pt() < params[7])
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
template <typename T>
bool selectionPIDKaon(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut)
{
  if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmatof_cut) {
    return true;
  }
  if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
    return true;
  }
  if (!use_tof && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
    return true;
  }
  return false;
}
template <typename T>
bool selectionPIDPion(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut)
{
  if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmatof_cut) {
    return true;
  }
  if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
    return true;
  }
  if (!use_tof && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
    return true;
  }
  return false;
}
template <typename T>
bool selectionPIDElec(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut)
{
  if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaEl() * candidate.tofNSigmaEl() + candidate.tpcNSigmaEl() * candidate.tpcNSigmaEl()) < nsigmatof_cut) {
    return true;
  }
  if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaEl()) < nsigmatpc_cut) {
    return true;
  }
  if (!use_tof && std::abs(candidate.tpcNSigmaEl()) < nsigmatpc_cut) {
    return true;
  }
  return false;
}
template <typename T>
bool selectionPIDProton(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut)
{
  if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < nsigmatof_cut) {
    return true;
  }
  if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmatpc_cut) {
    return true;
  }
  if (!use_tof && std::abs(candidate.tpcNSigmaPr()) < nsigmatpc_cut) {
    return true;
  }
  return false;
}
template <typename T>
bool selectionPIDMuon(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut)
{
  if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaMu() * candidate.tofNSigmaMu() + candidate.tpcNSigmaMu() * candidate.tpcNSigmaMu()) < nsigmatof_cut) {
    return true;
  }
  if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaMu()) < nsigmatpc_cut) {
    return true;
  }
  if (!use_tof && std::abs(candidate.tpcNSigmaMu()) < nsigmatpc_cut) {
    return true;
  }
  return false;
}
#endif // PWGUD_CORE_SGTRACKSELECTOR_H_
