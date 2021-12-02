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

#ifndef GFW__H
#define GFW__H
#include "GFWCumulant.h"
#include <vector>
#include <utility>
#include <algorithm>
#include "TString.h"
#include "TObjArray.h"
using std::vector;

class GFW
{
 public:
  struct Region {
    int Nhar, Npar, NpT;
    vector<int> NparVec;
    double EtaMin = -999;
    double EtaMax = -999;
    int BitMask = 1;
    TString rName = "";
    bool operator<(const Region& a) const
    {
      return EtaMin < a.EtaMin;
    };
    Region operator=(const Region& a)
    {
      Nhar = a.Nhar;
      Npar = a.Npar;
      NparVec = a.NparVec;
      NpT = a.NpT;
      EtaMin = a.EtaMin;
      EtaMax = a.EtaMax;
      rName = a.rName;
      BitMask = a.BitMask;
      return *this;
    };
    void PrintStructure() { printf("%s: eta [%f.. %f].", rName.Data(), EtaMin, EtaMax); };
  };
  struct CorrConfig {
    vector<vector<int>> Regs{};
    vector<vector<int>> Hars{};
    vector<int> Overlap;
    /*vector<int> Regs {};
    vector<int> Hars {};
    vector<int> Regs2 {};
    vector<int> Hars2 {};
    int Overlap1=-1;
    int Overlap2=-1;*/
    bool pTDif = kFALSE;
    TString Head = "";
  };
  GFW();
  ~GFW();
  vector<Region> fRegions;
  vector<GFWCumulant> fCumulants;
  vector<int> fEmptyInt;
  void AddRegion(TString refName, int lNhar, int lNpar, double lEtaMin, double lEtaMax, int lNpT = 1, int BitMask = 1);
  void AddRegion(TString refName, int lNhar, int* lNparVec, double lEtaMin, double lEtaMax, int lNpT = 1, int BitMask = 1);
  int CreateRegions();
  void Fill(double eta, int ptin, double phi, double weight, int mask, double secondWeight = -1);
  void Clear(); // { for(auto ptr = fCumulants.begin(); ptr!=fCumulants.end(); ++ptr) ptr->ResetQs(); };
  GFWCumulant GetCumulant(int index) { return fCumulants.at(index); };
  TComplex Calculate(TString config, bool SetHarmsToZero = kFALSE);
  CorrConfig GetCorrelatorConfig(TString config, TString head = "", bool ptdif = kFALSE);
  TComplex Calculate(CorrConfig corconf, int ptbin, bool SetHarmsToZero, bool DisableOverlap = kFALSE);

 private:
  bool fInitialized;
  void SplitRegions();
  GFWCumulant fEmptyCumulant;
  TComplex TwoRec(int n1, int n2, int p1, int p2, int ptbin, GFWCumulant*, GFWCumulant*, GFWCumulant*);
  TComplex RecursiveCorr(GFWCumulant* qpoi, GFWCumulant* qref, GFWCumulant* qol, int ptbin, vector<int>& hars, vector<int>& pows); // POI, Ref. flow, overlapping region
  TComplex RecursiveCorr(GFWCumulant* qpoi, GFWCumulant* qref, GFWCumulant* qol, int ptbin, vector<int>& hars);                    // POI, Ref. flow, overlapping region
  // Deprecated and not used (for now):
  void AddRegion(Region inreg) { fRegions.push_back(inreg); };
  Region GetRegion(int index) { return fRegions.at(index); };
  int FindRegionByName(TString refName);
  vector<TString> fCalculatedNames;
  vector<TComplex> fCalculatedQs;
  int FindCalculated(TString identifier);
  // Calculateing functions:
  TComplex Calculate(int poi, int ref, vector<int> hars, int ptbin = 0); // For differential, need POI and reference
  TComplex Calculate(int poi, vector<int> hars);                         // For integrated case
  // Process one string (= one region)
  TComplex CalculateSingle(TString config);

  bool SetHarmonicsToZero(TString& instr);
};
#endif