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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_STATHANDLER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_STATHANDLER_H_

#include "BootstrapProfile.h"

#include "TH1.h"
#include "TProfile.h"

#include <iostream>
#include <map>
#include <vector>

class TRandom;

class StatHandler
{
 public:
  StatHandler();
  explicit StatHandler(int n);
  virtual ~StatHandler();
  enum StatMethod {
    kBootstrap,
    kVarianceOfSamples,
    kJackknife
  };
  using StatFunc = TProfile* (StatHandler::*)(const std::vector<TH1D*>&);

  TH1D* getNominal(const std::vector<TH1D*>& samplehists, StatMethod method = kBootstrap);

 private:
  std::map<StatMethod, StatFunc> methodMap; //!
  int nSamples;
  int nxbins;
  std::vector<double> xbins;
  void setCurrentAxis(TH1D* tmph);
  TProfile* varianceOfSamples(const std::vector<TH1D*>& samplehists);
  TProfile* bootstrap(const std::vector<TH1D*>& samplehists);
  TProfile* jackknife(const std::vector<TH1D*>& samplehists);
  void applyStatErrors(TH1D* nominal, TProfile*& stat_profile);

  ClassDef(StatHandler, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_STATHANDLER_H_
