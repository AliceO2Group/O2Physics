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

/// \file FlowPtContainer.h
/// \brief Class to handle angular and transverse momentum correlations
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_

#include "GFW.h"
#include "GFWConfig.h"

#include <Framework/HistogramSpec.h>

#include <TCollection.h>
#include <TH1.h>
#include <TList.h>
#include <TNamed.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <variant>
#include <vector>

namespace o2::analysis::genericframework::eventweight
{
enum EventWeight {
  UnityWeight,
  TupleWeight
};
};

namespace o2::analysis::genericframework::flowptcontainer
{
struct PtProfileConfig {
  GFW::CorrConfig config;
  std::vector<std::vector<bool>> isPt;
  std::vector<std::vector<int>> factors;
  std::uint64_t activeFactors{}; // Flow, when present, is bit 0; pT terms follow in expression order.
  bool usesFullTopology{};
};

struct PtCorrConfig {
  std::string head;
  std::uint64_t allFactors{};
  bool hasFlow{};
  std::vector<PtProfileConfig> profiles;
  std::vector<std::string> meanProfiles;
};

struct PtRegionMoments {
  int maxOrder{};
  std::vector<double> values;
};

struct PtMomentGroup {
  std::string region;
  int order{};
  std::string meanProfile;
  ClassDefNV(PtMomentGroup, 1);
};

struct PtMomentConfig {
  std::string head;
  std::vector<PtMomentGroup> groups;
  ClassDefNV(PtMomentConfig, 1);
};
} // namespace o2::analysis::genericframework::flowptcontainer

class FlowPtContainer : public TNamed
{
 public:
  using FillType = std::variant<std::complex<double>, double>;
  FlowPtContainer();
  explicit FlowPtContainer(const char* name);
  ~FlowPtContainer() override;
  FlowPtContainer(const char* name, const char* title);
  void initialise(const o2::framework::AxisSpec& axis, const o2::analysis::genericframework::GFWRegions& regions, const o2::analysis::genericframework::GFWCorrConfigs& configs, const int& nsub = 10);
  void initialise(const o2::framework::AxisSpec& axis, const int& maxOrder, const o2::analysis::genericframework::GFWCorrConfigs& configs, const int& nsub = 10);
  void initialise(int nbinsx, double* xbins, const int& maxOrder, const o2::analysis::genericframework::GFWCorrConfigs& configs, const int& nsub = 10);
  void initialise(int nbinsx, double xlow, double xhigh, const int& maxOrder, const o2::analysis::genericframework::GFWCorrConfigs& configs, const int& nsub = 10);
  // initial pt-pt correlations with two subevents
  void initialiseSubevent(const o2::framework::AxisSpec& axis, const int& maxOrder, const int& nsubev = 2, const int& nsub = 10);
  void initialiseSubevent(int nbinsx, double* xbins, const int& maxOrder, const int& nsubev = 2, const int& nsub = 10);
  void initialiseSubevent(int nbinsx, double xlow, double xhigh, const int& maxOrder, const int& nsubev = 2, const int& nsub = 10);
  void fill(const double& w, const double& pt);
  void fill(const double& eta, const double& w, const double& pt, int mask);
  void fillSub(const double& w, const double& pt, int subIndex);
  void fillSub1(const double& w, const double& pt) { fillSub(w, pt, 0); }
  void fillSub2(const double& w, const double& pt) { fillSub(w, pt, nSubevents - 1); }
  void fillArray(FillType a, FillType b, double c, double d);
  int getVectorIndex(const int i, const int j) { return j * (mpar + 1) + i; }                                              // index for 2d array for storing pt correlations
  int getVectorIndex(const int i, const int j, const int k, const int l) { return i + j * 3 + k * 3 * 3 + l * 3 * 3 * 5; } // index for 4d array for std vnpt correlation - size 3x3x3x3
  void calculateCorrelations();
  void calculateSubeventCorrelations();
  void calculateCMTerms();
  void fillPtProfiles(const double& centmult, const double& rn);
  void fillPtProfiles(GFW& gfw, const double& centmult, const double& rn);
  // Book and fill a separate pT observable with a weight calculated by the task.
  // The task can supply a separate event weight for the same observable.
  bool addPtProfile(const char* name, int observableOrder);
  bool fillPtProfile(const char* name, int observableOrder, double mult, double eventWeight, double rn);
  void fillSubeventPtProfiles(const double& centmult, const double& rn);
  void fillVnPtCorrProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnDeltaPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnPtCorrProfiles(const int configIndex, const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnDeltaPtProfiles(const int configIndex, const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask);
  void fillVnDeltaPtStdProfiles(const double& centmult, const double& rn);
  void fillVnPtCorrStdProfiles(const double& centmult, const double& rn);
  void fillVnPtProfiles(const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
  {
    if (fUseCentralMoments) {
      fillVnDeltaPtProfiles(centmult, flowval, flowtuples, rn, mask);
    } else {
      fillVnPtCorrProfiles(centmult, flowval, flowtuples, rn, mask);
    }
  }
  void fillVnPtProfiles(const int configIndex, const double& centmult, const double& flowval, const double& flowtuples, const double& rn, uint8_t mask)
  {
    if (fUseCentralMoments) {
      fillVnDeltaPtProfiles(configIndex, centmult, flowval, flowtuples, rn, mask);
    } else {
      fillVnPtCorrProfiles(configIndex, centmult, flowval, flowtuples, rn, mask);
    }
  }
  void skipVnPtProfiles(uint8_t mask)
  {
    for (auto m(1); m <= mpar; ++m) {
      if ((mask & (1 << (m - 1))) == 0) {
        continue;
      }
      if (fUseCentralMoments) {
        for (auto i = 0; i <= m; ++i) {
          ++fillCounter;
        }
      } else {
        ++fillCounter;
      }
    }
  }
  void fillVnPtStdProfiles(const double& centmult, const double& rn)
  {
    if (fUseCentralMoments) {
      fillVnDeltaPtStdProfiles(centmult, rn);
    } else {
      fillVnPtCorrStdProfiles(centmult, rn);
    }
  }
  void fillCMProfiles(const double& centmult, const double& rn);
  void fillCMSubeventProfiles(const double& centmult, const double& rn);
  TList* getCorrList() { return fCorrList; }
  TList* getCMTermList() { return fCMTermList; }
  TList* getCovList() { return fCovList; }
  void setEventWeight(const unsigned int& lWeight) { fEventWeight = lWeight; }
  void setUseCentralMoments(bool newval) { fUseCentralMoments = newval; }
  void setUseGapMethod(bool newval) { fUseGap = newval; }
  bool usesCentralMoments() { return fUseCentralMoments; }
  bool usesGap() { return fUseGap; }
  void rebinMulti(int nbins);
  void rebinMulti(int nbins, double* binedges);
  TH1* getCentralMomentHist(int ind, int m);
  TH1* getCumulantHist(int ind, int m);
  TH1* getCorrHist(int ind, int m);
  // Postprocess a pure-pT configuration by head. ind = -1 selects the main profile.
  TH1* getPtCentralMomentHist(const char* head, int ind = -1);
  TH1* getPtCumulantHist(const char* head, int ind = -1);
  int getMpar() { return mpar; }
  Long64_t Merge(TCollection* collist);
  double orderedAddition(std::vector<double> vec);
  void createCentralMomentList();
  void calculateCentralMomentHists(std::vector<TH1*> inh, int ind, int m, TH1* hMpt);
  void createCumulantList();
  void calculateCumulantHists(std::vector<TH1*> inh, int ind);
  void clearVector()
  {
    sumP.clear();
    sumP.resize((mpar + 1) * (mpar + 1));
    insub.clear();
    insub.resize(nSubevents, std::vector<double>((mpar + 1) * (mpar + 1)));
    cmVal.clear();
    cmValSub.clear();
    cmValSub.resize(nSubevents);
    cmDen.clear();
    cmDenSub.clear();
    cmDenSub.resize(nSubevents);
    fillCounter = 0;
    arr.clear();
    arr.resize(3 * 3 * 5 * 5, {0.0, 0.0});
    warr.clear();
    warr.resize(3 * 3 * 5 * 5, 0.0);
    clearPtRegionMoments();
  };

  TList* fCMTermList;
  TList* fCorrList;
  TList* fCovList;
  TList* fSubList;
  TList* fSubCMList;
  TList* fCumulantList;
  TList* fCentralMomentList;

  int mpar;
  int nSubevents;            //!
  int fillCounter;           //!
  unsigned int fEventWeight; //!
  bool fUseCentralMoments;
  bool fUseGap;
  void mergeBSLists(TList* source, TList* target);
  TH1* raiseHistToPower(TH1* inh, double p);
  std::vector<double> sumP;                    //!
  std::vector<std::vector<double>> insub;      //!
  std::vector<double> corrNum;                 //!
  std::vector<std::vector<double>> corrNumSub; //!
  std::vector<double> corrDen;                 //!
  std::vector<std::vector<double>> corrDenSub; //!
  std::vector<double> cmVal;                   //!
  std::vector<std::vector<double>> cmValSub;   //!
  std::vector<double> cmDen;                   //!
  std::vector<std::vector<double>> cmDenSub;   //!
  std::vector<std::complex<double>> arr;       //!
  std::vector<double> warr;                    //!
  std::vector<int> fCovFirstIndex;             //!
  template <typename T>
  double getStdAABBCC(T& inarr);
  template <typename T>
  double getStdAABBCD(T& inarr);
  template <typename T>
  double getStdAABBDD(T& inarr);
  template <typename T>
  double getStdAABBC(T& inarr);
  template <typename T>
  double getStdAABBD(T& inarr);
  template <typename T>
  double getStdABCC(T& inarr);
  template <typename T>
  double getStdABCD(T& inarr);
  template <typename T>
  double getStdABDD(T& inarr);
  template <typename T>
  double getStdABC(T& inarr);
  template <typename T>
  double getStdABD(T& inarr);
  template <typename T>
  double getStdABCCCC(T& inarr);
  template <typename T>
  double getStdABCCCD(T& inarr);
  template <typename T>
  double getStdABCCDD(T& inarr);
  template <typename T>
  double getStdABCDDD(T& inarr);
  template <typename T>
  double getStdABDDDD(T& inarr);
  template <typename T>
  double getStdABCCC(T& inarr);
  template <typename T>
  double getStdABCCD(T& inarr);
  template <typename T>
  double getStdABCDD(T& inarr);
  template <typename T>
  double getStdABDDD(T& inarr);

 private:
  static const int centralMomentMaxOrder = 4;
  o2::analysis::genericframework::GFWRegions fPtRegions;                                          //!
  std::vector<o2::analysis::genericframework::flowptcontainer::PtCorrConfig> fPtCorrConfigs;      //!
  std::vector<o2::analysis::genericframework::flowptcontainer::PtRegionMoments> fPtRegionMoments; //!
  std::vector<o2::analysis::genericframework::flowptcontainer::PtMomentConfig> fPtMomentConfigs;
  TList* fPtCentralMomentList;             //!
  TList* fPtCumulantList;                  //!
  std::vector<std::vector<int>> subevents; //!
  void clearPtPostprocessLists();
  void clearPtRegionMoments();
  [[nodiscard]] double getPtMoment(int region, int weightPower, int ptPower) const;
  double recursivePtCorrelation(int poi, int ref, int overlap, std::vector<int>& weightPowers, std::vector<int>& ptPowers) const;
  void getSubevents(int k, int n, std::vector<int>& current, std::vector<std::vector<int>>& subevents);
  static constexpr std::array<float, 9> FactorialArray = {1., 1., 2., 6., 24., 120., 720., 5040., 40320.};
  static constexpr std::array<int, 9> SignArray = {1, -1, 1, -1, 1, -1, 1, -1, 1};
  ClassDefOverride(FlowPtContainer, 3);
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_FLOWPTCONTAINER_H_
