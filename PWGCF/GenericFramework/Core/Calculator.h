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

// NOLINTBEGIN(readability/braces)
#ifndef PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_
#define PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_

#include "BootstrapProfile.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"

#include <TH1.h>

#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

class Calculator
{
 public:
  using ObsFunc = TH1D* (Calculator::*)(int);
  using PtDiffFunc = TH1D* (Calculator::*)(int, double, double);
  Calculator();
  ~Calculator();
  bool readData(TFile* infile, std::string dirname, bool isGen = false, bool isRadialFlow = false);
  TH1D* calculate(ObsFunc func, int i);
  TH1D* calculate(PtDiffFunc func, int i, double min, double max);
  TH1D* calculateWithLateCentralityRebin(PtDiffFunc func, int i, double min, double max);
  void Rebin(int nrb, double* chbins);
  void setFlowID(std::string id) { flowid = id; }
  void setSystem(std::string system) { collision_system = system; }
  void setInvMassFitPath(std::string path)
  {
    invMassFitPath = path;
    invMassFitCache.clear();
  }
  void setInvMassFitFunctions(const std::string& species, const std::string& signalFunction, const std::string& backgroundFunction);
  void setInvMassFitFunctionsForSystematic(int systematicIndex,
                                           const std::string& species,
                                           const std::string& signalFunction,
                                           const std::string& backgroundFunction);
  void setInvMassFitCentralityBins(const std::vector<double>& centralityBins);
  void setUseTotal(bool newval) { useTotal = newval; }
  void setUseInverseVarianceWeight(bool newval) { useInverseVarianceWeight = newval; }

  struct RadialFlowAnalysis {
    std::vector<TH2D*> fraction_cent{4};        // ch, pi, ka, pr
    std::vector<TH2D*> fraction_cent_K0{3};     // sb1, sig, sb2
    std::vector<TH2D*> fraction_cent_Lambda{3}; // sb1, sig, sb2
  };

  // Flow functions
  TH1D* getVNN(int nhar, int npar, std::string id, int i);
  TH1D* getV12(int i) { return flowid.empty() ? getVNN(1, 2, "RefGap", i) : getVNN(1, 2, flowid, i); }
  TH1D* getV22(int i) { return flowid.empty() ? getVNN(2, 2, "RefGap", i) : getVNN(2, 2, flowid, i); }
  TH1D* getV32(int i) { return flowid.empty() ? getVNN(3, 2, "RefGap", i) : getVNN(3, 2, flowid, i); }
  TH1D* getV42(int i) { return flowid.empty() ? getVNN(4, 2, "RefGap", i) : getVNN(4, 2, flowid, i); }
  TH1D* getV24(int i) { return flowid.empty() ? getVNN(2, 4, "RefFull", i) : getVNN(2, 4, flowid, i); }
  TH1D* getV26(int i) { return flowid.empty() ? getVNN(2, 6, "RefFull", i) : getVNN(2, 6, flowid, i); }
  TH1D* getV28(int i) { return flowid.empty() ? getVNN(2, 8, "RefFull", i) : getVNN(2, 8, flowid, i); }
  TH1D* getCNN(int nhar, int npar, std::string id, int i);
  TH1D* getC12(int i) { return flowid.empty() ? getCNN(1, 2, "RefGap", i) : getCNN(1, 2, flowid, i); }
  TH1D* getC12Full(int i) { return getCNN(1, 2, "ChFull", i); }
  TH1D* getC22(int i) { return flowid.empty() ? getCNN(2, 2, "RefGap", i) : getCNN(2, 2, flowid, i); }
  TH1D* getC22Full(int i) { return getCNN(2, 2, "ChFull", i); }
  TH1D* getC32(int i) { return flowid.empty() ? getCNN(3, 2, "RefGap", i) : getCNN(3, 2, flowid, i); }
  TH1D* getC42(int i) { return flowid.empty() ? getCNN(4, 2, "RefGap", i) : getCNN(4, 2, flowid, i); }
  TH1D* getC24(int i) { return flowid.empty() ? getCNN(2, 4, "RefFull", i) : getCNN(2, 4, flowid, i); }
  TH1D* getC14(int i) { return flowid.empty() ? getCNN(1, 4, "RefGap", i) : getCNN(1, 4, flowid, i); }
  TH1D* getC26(int i) { return flowid.empty() ? getCNN(2, 6, "RefFull", i) : getCNN(2, 6, flowid, i); }
  TH1D* getC28(int i) { return flowid.empty() ? getCNN(2, 8, "RefFull", i) : getCNN(2, 8, flowid, i); }
  TH1D* getNN(int nhar, int npar, std::string id, int i);
  TH1D* get22(int i) { return flowid.empty() ? getNN(2, 2, "RefGap", i) : getNN(2, 2, flowid, i); }
  TH1D* get24(int i) { return flowid.empty() ? getNN(2, 4, "RefGap", i) : getNN(2, 4, flowid, i); }
  TH1D* get32(int i) { return flowid.empty() ? getNN(3, 2, "RefGap", i) : getNN(3, 2, flowid, i); }
  TH1D* get42(int i) { return flowid.empty() ? getNN(4, 2, "RefFull", i) : getNN(4, 2, flowid, i); }

  TH1D* getV32overV22(int i);
  TH1D* getV24overV22(int i);

  TH1D* getV24_3sub(int i);
  TH1D* getC24_3sub(int i);
  TH1D* getV2LLMR(int i);
  TH1D* getV2LMMR(int i);
  TH1D* getV2LMRR(int i);
  TH1D* getC2LLMR(int i);
  TH1D* getC2LMMR(int i);
  TH1D* getC2LMRR(int i);

  TH1D* get2LR(int i);
  TH1D* get2LM(int i);
  TH1D* get2MR(int i);

  // pt-diff flow
  TH1D* getVNN(int nhar, int npar, std::string id, int i, double min, double max);
  TH1D* getVPt22(int i, double min, double max) { return flowid.empty() ? getVNN(2, 2, "ChPGap", i, min, max) : getVNN(2, 2, flowid, i, min, max); }
  TH1D* getVPt32(int i, double min, double max) { return flowid.empty() ? getVNN(3, 2, "ChPGap", i, min, max) : getVNN(3, 2, flowid, i, min, max); }
  TH1D* getVPt42(int i, double min, double max) { return flowid.empty() ? getVNN(4, 2, "ChPGap", i, min, max) : getVNN(4, 2, flowid, i, min, max); }
  TH1D* getVPt24(int i, double min, double max) { return flowid.empty() ? getVNN(2, 4, "ChFull", i, min, max) : getVNN(2, 4, flowid, i, min, max); }
  TH1D* getCNN(int nhar, int npar, std::string id, int i, double min, double max);
  TH1D* getD22(int i, double min, double max) { return flowid.empty() ? getCNN(2, 2, "ChPGap", i, min, max) : getCNN(2, 2, flowid, i, min, max); }
  TH1D* getD32(int i, double min, double max) { return flowid.empty() ? getCNN(3, 2, "ChPGap", i, min, max) : getCNN(3, 2, flowid, i, min, max); }
  TH1D* getD42(int i, double min, double max) { return flowid.empty() ? getCNN(4, 2, "ChPGap", i, min, max) : getCNN(4, 2, flowid, i, min, max); }
  TH1D* getD24(int i, double min, double max) { return flowid.empty() ? getCNN(2, 4, "ChFull", i, min, max) : getCNN(2, 4, flowid, i, min, max); }

  TH1D* getVPt24_3sub(int i, double min, double max);
  TH1D* getVPt24_3sub_avgOnVn(int i, double min, double max);
  TH1D* getD24_3sub(int i, double min, double max);
  TH1D* getD2LpLMR(int i, double min, double max);
  TH1D* getD2LLMpR(int i, double min, double max);
  TH1D* getD2LLMRp(int i, double min, double max);
  TH1D* getD2LpMMR(int i, double min, double max);
  TH1D* getD2LMpMR(int i, double min, double max);
  TH1D* getD2LMMRp(int i, double min, double max);
  TH1D* getD2LpMRR(int i, double min, double max);
  TH1D* getD2LMpRR(int i, double min, double max);
  TH1D* getD2LMRpR(int i, double min, double max);

  TH1D* get2LpLMR(int i, double min, double max);
  TH1D* get2LLMpR(int i, double min, double max);
  TH1D* get2LLMRp(int i, double min, double max);
  TH1D* get2LpMMR(int i, double min, double max);
  TH1D* get2LMpMR(int i, double min, double max);
  TH1D* get2LMMRp(int i, double min, double max);
  TH1D* get2LpMRR(int i, double min, double max);
  TH1D* get2LMpRR(int i, double min, double max);
  TH1D* get2LMRpR(int i, double min, double max);
  TH1D* get2LRp(int i, double min, double max);
  TH1D* get2LpR(int i, double min, double max);
  TH1D* get2LpM(int i, double min, double max);
  TH1D* get2LMp(int i, double min, double max);
  TH1D* get2MpR(int i, double min, double max);
  TH1D* get2MRp(int i, double min, double max);

  // Flow variance
  TH1D* getVarV22(int i);

  // Flow-pt cumulants
  TH1D* getCVN2pt(int nhar, int i);
  TH1D* getCV12pt(int i) { return getCVN2pt(1, i); }
  TH1D* getCV22pt(int i) { return getCVN2pt(2, i); }
  TH1D* getCV32pt(int i) { return getCVN2pt(3, i); }
  TH1D* getCV42pt(int i) { return getCVN2pt(4, i); }
  // Could in theory be generalised to CVN2ptX, but will we ever need those? In any case trivial to change
  TH1D* getCV22pt2(int i);
  TH1D* getCV22pt3(int i);
  TH1D* getCV22pt4(int i);
  TH1D* getCV24pt(int i);
  TH1D* getCV24pt2(int i);

  // Flow-pt normalised cumulants
  TH1D* getRhoVN2pt(int nhar, int i);
  TH1D* getRhoV22pt(int i) { return getRhoVN2pt(2, i); }
  TH1D* getRhoV32pt(int i) { return getRhoVN2pt(3, i); }
  TH1D* getRhoV42pt(int i) { return getRhoVN2pt(4, i); }

  // Pt cumulant
  TH1D* getMpt(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCorrHist(i, 1) : flowptcont->getCorrHist(i, 1)); }
  TH1D* getMpt2(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCorrHist(i, 2) : flowptcont->getCorrHist(i, 2)); }
  TH1D* getMpt3(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCorrHist(i, 3) : flowptcont->getCorrHist(i, 3)); }
  TH1D* getMpt4(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCorrHist(i, 4) : flowptcont->getCorrHist(i, 4)); }

  TH1D* getK2(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCumulantHist(i, 2) : flowptcont->getCumulantHist(i, 2)); }
  TH1D* getK3(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCumulantHist(i, 3) : flowptcont->getCumulantHist(i, 3)); }
  TH1D* getK4(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCumulantHist(i, 4) : flowptcont->getCumulantHist(i, 4)); }
  TH1D* getMu2(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCentralMomentHist(i, 2) : flowptcont->getCentralMomentHist(i, 2)); }
  TH1D* getMu3(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCentralMomentHist(i, 3) : flowptcont->getCentralMomentHist(i, 3)); }
  TH1D* getMu4(int i) { return dynamic_cast<TH1D*>((flowptcontfull) ? flowptcontfull->getCentralMomentHist(i, 4) : flowptcont->getCentralMomentHist(i, 4)); }

  // Normalised pt cumulants
  TH1D* getK2Norm(int i);
  TH1D* getK3Norm(int i);
  TH1D* getDeltaPtOverPt(int i);
  TH1D* getStdSkew(int i);
  TH1D* getIntSkew(int i);
  TH1D* getExKur(int i);

  TH1D* getMptSubAvg(int i);
  // 2-subevent pt cumulants
  TH1D* getK2sub(int i);
  TH1D* getMu2sub(int i);
  TH1D* getK3sub(int i);
  TH1D* getK3subA(int i);
  TH1D* getK3subB(int i);
  TH1D* getMu3sub(int i);
  TH1D* getK4sub(int i);
  TH1D* getMu4sub(int i);

  TH1D* getK3_3sub(int i);
  TH1D* getK4_3sub(int i);

  // 4-subevent pt cumulants
  TH1D* getK4_4sub(int i);

  // normalised subevent pt cumulants
  TH1D* getK2NormSub(int i);
  TH1D* getK3NormSub(int i);
  TH1D* getStdSkewSub(int i);
  TH1D* getIntSkewSub(int i);
  TH1D* getIntSkewAvgSub(int i);
  TH1D* getExKurSub(int i);
  TH1D* getK3Norm_3sub(int i);
  TH1D* getStdSkew_3sub(int i);
  TH1D* getIntSkew_3sub(int i);
  TH1D* getExKur_4sub(int i);

  TH1D* getIntegratedV0(int i);
  TH1D* getIntegratedV02(int i);
  // Radial flow analysis
  TH1D* getnpt(int i, int pidIndex, double min, double max);
  TH1D* getnptCh(int i, double min, double max) { return getnpt(i, 0, min, max); };
  TH1D* getnptPi(int i, double min, double max) { return getnpt(i, 1, min, max); };
  TH1D* getnptKa(int i, double min, double max) { return getnpt(i, 2, min, max); };
  TH1D* getnptPr(int i, double min, double max) { return getnpt(i, 3, min, max); };
  TH1D* getnptK0(int i, double min, double max);
  TH1D* getnptLambda(int i, double min, double max);
  TH1D* getV0(int i, int pidIndex, double min, double max);
  TH1D* getV0Ch(int i, double min, double max) { return getV0(i, 0, min, max); };
  TH1D* getV0Pi(int i, double min, double max) { return getV0(i, 1, min, max); };
  TH1D* getV0Ka(int i, double min, double max) { return getV0(i, 2, min, max); };
  TH1D* getV0Pr(int i, double min, double max) { return getV0(i, 3, min, max); };
  TH1D* getV0K0(int i, double min, double max);
  TH1D* getV0Lambda(int i, double min, double max);
  // Four-term centred numerators. The existing getV0/getV02 methods remain available.
  // pidIndex: 0=ch, 1=pi, 2=ka, 3=pr, 4=K0, 5=Lambda.
  // Requires the _w2pc/_w3pc fractions and ChGap3pcW22 from the weighted flow task.
  // With cfgUsePIDTotal=true, V02 for identified species and V0 for K0/Lambda
  // need additional species-specific combined-weight profiles.
  TH1D* getV0NumeratorWeighted(int i, int pidIndex, double min, double max);
  TH1D* getV02NumeratorWeighted(int i, int pidIndex, double min, double max);
  TH1D* getV0Weighted(int i, int pidIndex, double min, double max);
  TH1D* getV00Weighted(int i, int pidIndex, double min, double max);
  TH1D* getV02Weighted(int i, int pidIndex, double min, double max);
  TH1D* getV0WeightedCh(int i, double min, double max) { return getV0Weighted(i, 0, min, max); }
  TH1D* getV0WeightedPi(int i, double min, double max) { return getV0Weighted(i, 1, min, max); }
  TH1D* getV0WeightedKa(int i, double min, double max) { return getV0Weighted(i, 2, min, max); }
  TH1D* getV0WeightedPr(int i, double min, double max) { return getV0Weighted(i, 3, min, max); }
  TH1D* getV0WeightedK0(int i, double min, double max) { return getV0Weighted(i, 4, min, max); }
  TH1D* getV0WeightedLambda(int i, double min, double max) { return getV0Weighted(i, 5, min, max); }
  // New definition
  TH1D* getV00WeightedCh(int i, double min, double max) { return getV00Weighted(i, 0, min, max); }
  TH1D* getV00WeightedPi(int i, double min, double max) { return getV00Weighted(i, 1, min, max); }
  TH1D* getV00WeightedKa(int i, double min, double max) { return getV00Weighted(i, 2, min, max); }
  TH1D* getV00WeightedPr(int i, double min, double max) { return getV00Weighted(i, 3, min, max); }
  TH1D* getV00WeightedK0(int i, double min, double max) { return getV00Weighted(i, 4, min, max); }
  TH1D* getV00WeightedLambda(int i, double min, double max) { return getV00Weighted(i, 5, min, max); }
  TH1D* getV02WeightedCh(int i, double min, double max) { return getV02Weighted(i, 0, min, max); }
  TH1D* getV02WeightedPi(int i, double min, double max) { return getV02Weighted(i, 1, min, max); }
  TH1D* getV02WeightedKa(int i, double min, double max) { return getV02Weighted(i, 2, min, max); }
  TH1D* getV02WeightedPr(int i, double min, double max) { return getV02Weighted(i, 3, min, max); }
  TH1D* getV02WeightedK0(int i, double min, double max) { return getV02Weighted(i, 4, min, max); }
  TH1D* getV02WeightedLambda(int i, double min, double max) { return getV02Weighted(i, 5, min, max); }
  TH1D* getV02(int i, int pidIndex, double min, double max);
  TH1D* getV02Ch(int i, double min, double max) { return getV02(i, 0, min, max); };
  TH1D* getV02Pi(int i, double min, double max) { return getV02(i, 1, min, max); };
  TH1D* getV02Ka(int i, double min, double max) { return getV02(i, 2, min, max); };
  TH1D* getV02Pr(int i, double min, double max) { return getV02(i, 3, min, max); };
  TH1D* getV02K0(int i, double min, double max);
  TH1D* getV02Lambda(int i, double min, double max);
  TH1D* getV02K0Alt(int i, double min, double max);
  TH1D* getV02LambdaAlt(int i, double min, double max);
  TH1D* getV0K0Alt(int i, double min, double max);
  TH1D* getV0LambdaAlt(int i, double min, double max);

  TH1D* getV0Corr(int i, int pidIndex, double min, double max);
  TH1D* getV0CorrCh(int i, double min, double max) { return getV0Corr(i, 0, min, max); };
  TH1D* getV0CorrPi(int i, double min, double max) { return getV0Corr(i, 1, min, max); };
  TH1D* getV0CorrKa(int i, double min, double max) { return getV0Corr(i, 2, min, max); };
  TH1D* getV0CorrPr(int i, double min, double max) { return getV0Corr(i, 3, min, max); };
  TH1D* getV0CorrK0(int i, double min, double max);
  TH1D* getV0CorrLambda(int i, double min, double max);
  TH1D* getV0Cov(int i, int pidIndex, double min, double max);
  TH1D* getV0CovCh(int i, double min, double max) { return getV02Cov(i, 0, min, max); };
  TH1D* getV0CovPi(int i, double min, double max) { return getV02Cov(i, 1, min, max); };
  TH1D* getV0CovKa(int i, double min, double max) { return getV02Cov(i, 2, min, max); };
  TH1D* getV0CovPr(int i, double min, double max) { return getV02Cov(i, 3, min, max); };
  TH1D* getV0CovK0(int i, double min, double max);
  TH1D* getV0CovLambda(int i, double min, double max);
  TH1D* getV02Corr(int i, int pidIndex, double min, double max);
  TH1D* getV02CorrCh(int i, double min, double max) { return getV02Corr(i, 0, min, max); };
  TH1D* getV02CorrPi(int i, double min, double max) { return getV02Corr(i, 1, min, max); };
  TH1D* getV02CorrKa(int i, double min, double max) { return getV02Corr(i, 2, min, max); };
  TH1D* getV02CorrPr(int i, double min, double max) { return getV02Corr(i, 3, min, max); };
  TH1D* getV02CorrK0(int i, double min, double max);
  TH1D* getV02CorrLambda(int i, double min, double max);
  TH1D* getV02Cov(int i, int pidIndex, double min, double max);
  TH1D* getV02CovCh(int i, double min, double max) { return getV02Cov(i, 0, min, max); };
  TH1D* getV02CovPi(int i, double min, double max) { return getV02Cov(i, 1, min, max); };
  TH1D* getV02CovKa(int i, double min, double max) { return getV02Cov(i, 2, min, max); };
  TH1D* getV02CovPr(int i, double min, double max) { return getV02Cov(i, 3, min, max); };
  TH1D* getV02CovK0(int i, double min, double max);
  TH1D* getV02CovLambda(int i, double min, double max);

  TH1D* getV0Means(int i, int pidIndex, double min, double max);
  TH1D* getV0MeansCh(int i, double min, double max) { return getV0Means(i, 0, min, max); };
  TH1D* getV0MeansPi(int i, double min, double max) { return getV0Means(i, 1, min, max); };
  TH1D* getV0MeansKa(int i, double min, double max) { return getV0Means(i, 2, min, max); };
  TH1D* getV0MeansPr(int i, double min, double max) { return getV0Means(i, 3, min, max); };
  TH1D* getV0MeansK0(int i, double min, double max);
  TH1D* getV0MeansLambda(int i, double min, double max);
  TH1D* getV02Means(int i, int pidIndex, double min, double max);
  TH1D* getV02MeansCh(int i, double min, double max) { return getV02Means(i, 0, min, max); };
  TH1D* getV02MeansPi(int i, double min, double max) { return getV02Means(i, 1, min, max); };
  TH1D* getV02MeansKa(int i, double min, double max) { return getV02Means(i, 2, min, max); };
  TH1D* getV02MeansPr(int i, double min, double max) { return getV02Means(i, 3, min, max); };
  TH1D* getV02MeansK0(int i, double min, double max);
  TH1D* getV02MeansLambda(int i, double min, double max);

  // Two-component model
  TH1D* getAlpha(int i, double min, double max);
  TH1D* getBeta(int i, double min, double max);
  // Based on <[mpt]> normalised v00
  TH1D* getAlphaAlt(int i, double min, double max);
  TH1D* getBetaAlt(int i, double min, double max);

  TH1D* getCorrVsPt(int nhar, int npar, std::string id, int i, double min, double max);

  double getIntegratedValue(TH1D* inhist, double min, double max);
  double readInvMassFit(std::string species, double massMin = -1., double massMax = -1., double pTMin = 0.2, double pTMax = 3.0, int centMin = 0, int centMax = 20);

  TH1D* powerHistogram(TH1D* inhist, double power);

 private:
  enum class FractionScheme { V02,
                              V0 };
  // Data members
  FlowContainer* flowcont;                                                                                   //!
  FlowPtContainer* flowptcont;                                                                               //!
  FlowPtContainer* flowptcontfull;                                                                           //!
  RadialFlowAnalysis radialFlowV02;                                                                          //!
  RadialFlowAnalysis radialFlowV0;                                                                           //!
  RadialFlowAnalysis radialFlowV02Weighted;                                                                  //!
  RadialFlowAnalysis radialFlowV0Weighted;                                                                   //!
  TList* corr_list;                                                                                          //!
  std::string flowid;                                                                                        //!
  std::string collision_system;                                                                              //!
  std::string invMassFitPath;                                                                                //!
  std::map<std::string, std::pair<std::string, std::string>> invMassFitFunctions;                            //!
  std::map<int, std::map<std::string, std::pair<std::string, std::string>>> invMassFitFunctionsBySystematic; //!
  std::vector<double> invMassFitCentralityBins;                                                              //!
  int systematicIndex = 0;                                                                                   //!
  int invMassFitLookupCentMin = -1;                                                                          //!
  int invMassFitLookupCentMax = -1;                                                                          //!
  std::map<std::string, double> invMassFitCache;                                                             //!
  std::set<std::string> invMassFitFallbackWarnings;                                                          //!
  FlowContainer* activeFlowContainer = nullptr;                                                              //!
  int activeFlowContainerSample = -2;                                                                        //!
  bool hasChGap3pcW22 = false;                                                                               //!
  std::pair<int, int> subeventInfo;                                                                          //!
  bool is_digits(std::string& str);
  bool useTotal;
  bool useInverseVarianceWeight;
  std::string padWithZeros(int number, int totalDigits);
  std::pair<int, int> checkSubeventInput(FlowPtContainer* cont);
  void loadFractionHists(TFile* infile, std::string dirname, std::string suffix, RadialFlowAnalysis& radialFlowStruct, const std::string& profileSuffix = "");
  void clearFractionHists(RadialFlowAnalysis& radialFlowStruct);
  void clearInputData();
  std::unique_ptr<FlowContainer> prepareFlowContainerSample(int sampleIndex);
  FlowContainer* getFlowContainerSample(int sampleIndex, std::unique_ptr<FlowContainer>& ownedContainer);
  TH2D* getNptFile(TFile* infile, const std::string& dirname, const std::string& suffix, const std::string& species, const std::string& fallbackSuffix = "", const std::string& profileSuffix = "");
  TH1D* getProjectedFractionalSpectra(std::string sp, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getProjectedFractionalK0Spectra(std::string sp, int resoIndex, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getProjectedFractionalLambdaSpectra(std::string sp, int resoIndex, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getnpt(int i, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getnptK0(int i, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getnptLambda(int i, double min, double max, FractionScheme scheme, bool combinedWeight = false);
  TH1D* getFractionForSpecies(int i, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight);
  TH1D* getCorrelationForSpecies(int i, int pidIndex, double min, double max, FractionScheme scheme);
  ClassDef(Calculator, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_
// NOLINTEND(readability/braces)
