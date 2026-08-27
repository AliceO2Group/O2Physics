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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_SYSHANDLER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_SYSHANDLER_H_

#include "TH1.h"

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

class SysHandler
{
 public:
  SysHandler();
  explicit SysHandler(int n);
  virtual ~SysHandler();

  // Main function call
  TH1D* getNominal(bool relative = kTRUE, bool apply_barlow = kTRUE);

  // Utility functions
  struct SystSettings {
    std::pair<double, double> sys_avg_range;
    bool use_fit_for_syst;
    double max_allowed_syst;
    double fNbarlow;
    bool fUseBarlowCount;
  };
  struct SystematicVariation {
    int index;
    int group;
    int correlation;
    std::string description;
    std::string groupDescription;
  };
  void setSystSettings(SystSettings settings) { current_settings = settings; }
  const SystSettings getSystSettings() const { return current_settings; }
  void setSystAverageRange(std::pair<double, double> range, double lMaxAllowed = 1)
  {
    current_settings.sys_avg_range = range;
    current_settings.max_allowed_syst = lMaxAllowed;
  };
  void setBarlowCriteria(int newval) { current_settings.fNbarlow = newval; }
  void setUseBarlowCount(bool newval) { current_settings.fUseBarlowCount = newval; }
  void setSysMap(std::map<int, TH1*> sysmap)
  {
    clearSysMap();
    fSystMap = std::move(sysmap);
    fIndMap.clear();
  }
  void setSystematicVariations(const std::vector<SystematicVariation>& variations);
  std::map<int, double> getPrintMap() { return fPrintMap; }
  std::map<int, std::vector<double>> getPrintBinMap() { return fPrintBinMap; }
  void histOper(TH1* reth, Double_t cvl = 0, Double_t lpower = 1, Bool_t Abs = kFALSE, Bool_t RemoveErrors = kFALSE);
  std::map<int, TH1*> getSystMap() { return fSystMap; }
  TH1* getVariation(int index) { return fSystMap[index]; }

  // Systematics calculator functions
  TH1* getSystematics();
  TH1* getBarlowTest(int sys_index, double corrFactor = 1);
  TH1* getSystematics(int sys_index, bool relative = kTRUE, bool apply_barlow = kFALSE, bool abs = kTRUE);
  TH1* getSystematicsFromFit(int sys_index, bool apply_barlow = kFALSE, bool ratio = kFALSE, bool relative = kTRUE);
  TH1* getSystematicsSummed(bool relative = kTRUE, bool apply_barlow = kFALSE);
  std::vector<TH1*> getMergedErrors(bool relative = kTRUE, bool apply_barlow = kFALSE);
  std::vector<TH1*> getSystSubset(Int_t KeyInd, bool relative = kTRUE, bool apply_barlow = kFALSE);
  void applyErrors(TH1* inh, TH1* errors, Bool_t relative = kTRUE);
  void applyBarlow(TH1* target, TH1* barlow);
  void addInQuad(TH1* h1, TH1* h2);
  void RemoveErrors(TH1* inh)
  {
    for (Int_t i = 1; i <= inh->GetNbinsX(); i++)
      if (inh->GetBinError(i))
        inh->SetBinError(i, 0);
  };
  std::vector<TH1*> nullvec;
  int getCorrelationError(int fSystInd) const;
  std::string getSystDescription(int fSystInd) const;
  std::string getSystGroupDescription(int targetInd) const;
  int getTarget(int fSystInd) const;

 private:
  std::map<int, std::vector<int>> fIndMap;                  //! map to hold correlated indices
  std::map<int, SystematicVariation> fSystematicVariations; //! systematic metadata indexed by variation ID
  std::map<int, TH1*> fSystMap;                             //! map to hold variations
  std::map<int, double> fSystOverride;                      //! map to count correlated observables for error merging
  std::map<int, double> fPrintMap;                          //! map to store values for printing
  std::map<int, std::vector<double>> fPrintBinMap;          //! map to store bin-by-bin values for printing
  SystSettings current_settings;

  bool BuildIndexMap(Bool_t force = kFALSE);
  bool passesBarlowSelection(const TH1* barlow) const;
  void clearSysMap();

  ClassDef(SysHandler, 2);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_SYSHANDLER_H_
