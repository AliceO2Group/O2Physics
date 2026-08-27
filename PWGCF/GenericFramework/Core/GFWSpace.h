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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GFWSPACE_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GFWSPACE_H_

#include "Calculator.h"
namespace GFWSpace
{
Calculator::ObsFunc kNone = 0;
Calculator::PtDiffFunc kNonePt = 0;

// Flow observables
Calculator::ObsFunc kV12 = &Calculator::getV12; // v1{2}
Calculator::ObsFunc kV22 = &Calculator::getV22; // v2{2}
Calculator::ObsFunc kV32 = &Calculator::getV32; // v3{2}
Calculator::ObsFunc kV42 = &Calculator::getV42; // v4{2}

Calculator::ObsFunc kV24 = &Calculator::getV24;           // v2{4}
Calculator::ObsFunc kV26 = &Calculator::getV26;           // v2{6}
Calculator::ObsFunc kV28 = &Calculator::getV28;           // v2{8}
Calculator::ObsFunc kV24_3sub = &Calculator::getV24_3sub; // v2{4}_3sub
Calculator::ObsFunc kV2LLMR = &Calculator::getV2LLMR;     // v2{LLMR}
Calculator::ObsFunc kV2LMMR = &Calculator::getV2LMMR;     // v2{LMMR}
Calculator::ObsFunc kV2LMRR = &Calculator::getV2LMRR;     // v2{LMRR}

Calculator::ObsFunc kC12 = &Calculator::getC12;           // c1{2}
Calculator::ObsFunc kC12Full = &Calculator::getC12Full;   // c1{2}
Calculator::ObsFunc kC22 = &Calculator::getC22;           // c2{2}
Calculator::ObsFunc kC22Full = &Calculator::getC22Full;   // c2{2}
Calculator::ObsFunc kC32 = &Calculator::getC32;           // c3{2}
Calculator::ObsFunc kC42 = &Calculator::getC42;           // c4{2}
Calculator::ObsFunc kC24 = &Calculator::getC24;           // c2{4}
Calculator::ObsFunc kC14 = &Calculator::getC14;           // c1{4}
Calculator::ObsFunc kC26 = &Calculator::getC26;           // c2{6}
Calculator::ObsFunc kC28 = &Calculator::getC28;           // c2{8}
Calculator::ObsFunc kC24_3sub = &Calculator::getC24_3sub; // v2{4}_3sub
Calculator::ObsFunc kC2LLMR = &Calculator::getC2LLMR;     // c2{LLMR}
Calculator::ObsFunc kC2LMMR = &Calculator::getC2LMMR;     // c2{LMMR}
Calculator::ObsFunc kC2LMRR = &Calculator::getC2LMRR;     // c2{LMRR}

Calculator::ObsFunc k2LR = &Calculator::get2LR;
Calculator::ObsFunc k2LM = &Calculator::get2LM;
Calculator::ObsFunc k2MR = &Calculator::get2MR;

Calculator::ObsFunc k22 = &Calculator::get22; //<2>_{2,-2}
Calculator::ObsFunc k32 = &Calculator::get32; //<2>_{3,-3}
Calculator::ObsFunc k42 = &Calculator::get42; //<2>_{4,-4}
Calculator::ObsFunc k24 = &Calculator::get24; //<4>_{2,2,-2,-2}

Calculator::ObsFunc kV32overV22 = &Calculator::getV32overV22; // v3{2}/v2{2}
Calculator::ObsFunc kV24overV22 = &Calculator::getV24overV22; // v2{4}/v2{2}

Calculator::ObsFunc kIntegratedV0 = &Calculator::getIntegratedV0;
Calculator::ObsFunc kIntegratedV02 = &Calculator::getIntegratedV02;

// pt-diff flow
Calculator::PtDiffFunc kVPt22 = &Calculator::getVPt22; // pt-diff v2{2}
Calculator::PtDiffFunc kVPt32 = &Calculator::getVPt32; // pt-diff v3{2}
Calculator::PtDiffFunc kVPt42 = &Calculator::getVPt42; // pt-diff v4{2}
Calculator::PtDiffFunc kVPt24 = &Calculator::getVPt24; // pt-diff v2{4}
Calculator::PtDiffFunc kD22 = &Calculator::getD22;     // pt-diff v2{2}
Calculator::PtDiffFunc kD32 = &Calculator::getD32;     // pt-diff v3{2}
Calculator::PtDiffFunc kD42 = &Calculator::getD42;     // pt-diff v4{2}
Calculator::PtDiffFunc kD24 = &Calculator::getD24;     // pt-diff v2{4}

Calculator::PtDiffFunc kVPt24_3sub = &Calculator::getVPt24_3sub;
Calculator::PtDiffFunc kVPt24_3sub_avgvn = &Calculator::getVPt24_3sub_avgOnVn; // pt-diff v2{4}
Calculator::PtDiffFunc kD24_3sub = &Calculator::getD24_3sub;
Calculator::PtDiffFunc kD2LpLMR = &Calculator::getD2LpLMR;
Calculator::PtDiffFunc kD2LLMpR = &Calculator::getD2LLMpR;
Calculator::PtDiffFunc kD2LLMRp = &Calculator::getD2LLMRp;
Calculator::PtDiffFunc kD2LpMMR = &Calculator::getD2LpMMR;
Calculator::PtDiffFunc kD2LMpMR = &Calculator::getD2LMpMR;
Calculator::PtDiffFunc kD2LMMRp = &Calculator::getD2LMMRp;
Calculator::PtDiffFunc kD2LpMRR = &Calculator::getD2LpMRR;
Calculator::PtDiffFunc kD2LMpRR = &Calculator::getD2LMpRR;
Calculator::PtDiffFunc kD2LMRpR = &Calculator::getD2LMRpR;

Calculator::PtDiffFunc k2LpLMR = &Calculator::get2LpLMR;
Calculator::PtDiffFunc k2LLMpR = &Calculator::get2LLMpR;
Calculator::PtDiffFunc k2LLMRp = &Calculator::get2LLMRp;
Calculator::PtDiffFunc k2LpMMR = &Calculator::get2LpMMR;
Calculator::PtDiffFunc k2LMpMR = &Calculator::get2LMpMR;
Calculator::PtDiffFunc k2LMMRp = &Calculator::get2LMMRp;
Calculator::PtDiffFunc k2LpMRR = &Calculator::get2LpMRR;
Calculator::PtDiffFunc k2LMpRR = &Calculator::get2LMpRR;
Calculator::PtDiffFunc k2LMRpR = &Calculator::get2LMRpR;

Calculator::PtDiffFunc k2LRp = &Calculator::get2LRp;
Calculator::PtDiffFunc k2LpR = &Calculator::get2LpR;
Calculator::PtDiffFunc k2LpM = &Calculator::get2LpM;
Calculator::PtDiffFunc k2LMp = &Calculator::get2LMp;
Calculator::PtDiffFunc k2MpR = &Calculator::get2MpR;
Calculator::PtDiffFunc k2MRp = &Calculator::get2MRp;

Calculator::PtDiffFunc knptCh = &Calculator::getnptCh;
Calculator::PtDiffFunc knptPi = &Calculator::getnptPi;
Calculator::PtDiffFunc knptKa = &Calculator::getnptKa;
Calculator::PtDiffFunc knptPr = &Calculator::getnptPr;
Calculator::PtDiffFunc knptK0 = &Calculator::getnptK0;
Calculator::PtDiffFunc knptLambda = &Calculator::getnptLambda;

Calculator::PtDiffFunc kV0CorrCh = &Calculator::getV0CorrCh;
Calculator::PtDiffFunc kV0CorrPi = &Calculator::getV0CorrPi;
Calculator::PtDiffFunc kV0CorrKa = &Calculator::getV0CorrKa;
Calculator::PtDiffFunc kV0CorrPr = &Calculator::getV0CorrPr;
Calculator::PtDiffFunc kV0CorrK0 = &Calculator::getV0CorrK0;
Calculator::PtDiffFunc kV0CorrLambda = &Calculator::getV0CorrLambda;

Calculator::PtDiffFunc kV0CovCh = &Calculator::getV0CovCh;
Calculator::PtDiffFunc kV0CovPi = &Calculator::getV0CovPi;
Calculator::PtDiffFunc kV0CovKa = &Calculator::getV0CovKa;
Calculator::PtDiffFunc kV0CovPr = &Calculator::getV0CovPr;
Calculator::PtDiffFunc kV0CovK0 = &Calculator::getV0CovK0;
Calculator::PtDiffFunc kV0CovLambda = &Calculator::getV0CovLambda;

Calculator::PtDiffFunc kV0MeansCh = &Calculator::getV0MeansCh;
Calculator::PtDiffFunc kV0MeansPi = &Calculator::getV0MeansPi;
Calculator::PtDiffFunc kV0MeansKa = &Calculator::getV0MeansKa;
Calculator::PtDiffFunc kV0MeansPr = &Calculator::getV0MeansPr;
Calculator::PtDiffFunc kV0MeansK0 = &Calculator::getV0MeansK0;
Calculator::PtDiffFunc kV0MeansLambda = &Calculator::getV0MeansLambda;

Calculator::PtDiffFunc kV0Ch = &Calculator::getV0Ch;
Calculator::PtDiffFunc kV0Pi = &Calculator::getV0Pi;
Calculator::PtDiffFunc kV0Ka = &Calculator::getV0Ka;
Calculator::PtDiffFunc kV0Pr = &Calculator::getV0Pr;
Calculator::PtDiffFunc kV0K0 = &Calculator::getV0K0;
Calculator::PtDiffFunc kV0Lambda = &Calculator::getV0Lambda;
Calculator::PtDiffFunc kV0WeightedCh = &Calculator::getV0WeightedCh;
Calculator::PtDiffFunc kV0WeightedPi = &Calculator::getV0WeightedPi;
Calculator::PtDiffFunc kV0WeightedKa = &Calculator::getV0WeightedKa;
Calculator::PtDiffFunc kV0WeightedPr = &Calculator::getV0WeightedPr;
Calculator::PtDiffFunc kV0WeightedK0 = &Calculator::getV0WeightedK0;
Calculator::PtDiffFunc kV0WeightedLambda = &Calculator::getV0WeightedLambda;
Calculator::PtDiffFunc kV0K0Alt = &Calculator::getV0K0Alt;
Calculator::PtDiffFunc kV0LambdaAlt = &Calculator::getV0LambdaAlt;

Calculator::PtDiffFunc kV02CorrCh = &Calculator::getV02CorrCh;
Calculator::PtDiffFunc kV02CorrPi = &Calculator::getV02CorrPi;
Calculator::PtDiffFunc kV02CorrKa = &Calculator::getV02CorrKa;
Calculator::PtDiffFunc kV02CorrPr = &Calculator::getV02CorrPr;
Calculator::PtDiffFunc kV02CorrK0 = &Calculator::getV02CorrK0;
Calculator::PtDiffFunc kV02CorrLambda = &Calculator::getV02CorrLambda;

Calculator::PtDiffFunc kV02CovCh = &Calculator::getV02CovCh;
Calculator::PtDiffFunc kV02CovPi = &Calculator::getV02CovPi;
Calculator::PtDiffFunc kV02CovKa = &Calculator::getV02CovKa;
Calculator::PtDiffFunc kV02CovPr = &Calculator::getV02CovPr;
Calculator::PtDiffFunc kV02CovK0 = &Calculator::getV02CovK0;
Calculator::PtDiffFunc kV02CovLambda = &Calculator::getV02CovLambda;

Calculator::PtDiffFunc kV02MeansCh = &Calculator::getV02MeansCh;
Calculator::PtDiffFunc kV02MeansPi = &Calculator::getV02MeansPi;
Calculator::PtDiffFunc kV02MeansKa = &Calculator::getV02MeansKa;
Calculator::PtDiffFunc kV02MeansPr = &Calculator::getV02MeansPr;
Calculator::PtDiffFunc kV02MeansK0 = &Calculator::getV02MeansK0;
Calculator::PtDiffFunc kV02MeansLambda = &Calculator::getV02MeansLambda;

Calculator::PtDiffFunc kV02Ch = &Calculator::getV02Ch;
Calculator::PtDiffFunc kV02Pi = &Calculator::getV02Pi;
Calculator::PtDiffFunc kV02Ka = &Calculator::getV02Ka;
Calculator::PtDiffFunc kV02Pr = &Calculator::getV02Pr;
Calculator::PtDiffFunc kV02K0 = &Calculator::getV02K0;
Calculator::PtDiffFunc kV02Lambda = &Calculator::getV02Lambda;
Calculator::PtDiffFunc kV02WeightedCh = &Calculator::getV02WeightedCh;
Calculator::PtDiffFunc kV02WeightedPi = &Calculator::getV02WeightedPi;
Calculator::PtDiffFunc kV02WeightedKa = &Calculator::getV02WeightedKa;
Calculator::PtDiffFunc kV02WeightedPr = &Calculator::getV02WeightedPr;
Calculator::PtDiffFunc kV02WeightedK0 = &Calculator::getV02WeightedK0;
Calculator::PtDiffFunc kV02WeightedLambda = &Calculator::getV02WeightedLambda;
Calculator::PtDiffFunc kV02K0Alt = &Calculator::getV02K0Alt;
Calculator::PtDiffFunc kV02LambdaAlt = &Calculator::getV02LambdaAlt;

Calculator::ObsFunc kVarV22 = &Calculator::getVarV22; // Var(v_2^2)

// Transverse momentum observables
Calculator::ObsFunc kMpt = &Calculator::getMpt;   //[pt]
Calculator::ObsFunc kMpt2 = &Calculator::getMpt2; // 2nd moment of [pt]
Calculator::ObsFunc kMpt3 = &Calculator::getMpt3; // 3rd moment of [pt]
Calculator::ObsFunc kMpt4 = &Calculator::getMpt4; // 4th moment of [pt]
Calculator::ObsFunc kK2 = &Calculator::getK2;     // 2nd order [pt] cumulant
Calculator::ObsFunc kK3 = &Calculator::getK3;     // 3rd order [pt] cumulant
Calculator::ObsFunc kK4 = &Calculator::getK4;     // 4th order [pt] cumulant
Calculator::ObsFunc kMu2 = &Calculator::getMu2;   // 2nd order central [pt] moment (also known as ck)
Calculator::ObsFunc kMu3 = &Calculator::getMu3;   // 3rd order central [pt] moment
Calculator::ObsFunc kMu4 = &Calculator::getMu4;   // 4th order central [pt] moment

Calculator::ObsFunc kK2Norm = &Calculator::getK2Norm;               // normalised [pt] variance
Calculator::ObsFunc kDeltaPtOverPt = &Calculator::getDeltaPtOverPt; // normalised [pt] variance
Calculator::ObsFunc kK3Norm = &Calculator::getK3Norm;               // normalised [pt] variance
Calculator::ObsFunc kStdSkew = &Calculator::getStdSkew;             // standardized [pt] skewness
Calculator::ObsFunc kIntSkew = &Calculator::getIntSkew;             // intensive [pt] skewness
Calculator::ObsFunc kExKur = &Calculator::getExKur;                 // Excess [pt] kurtosis

Calculator::ObsFunc kMptSubAvg = &Calculator::getMptSubAvg;         //[pt]
Calculator::ObsFunc kK2Sub = &Calculator::getK2sub;                 //
Calculator::ObsFunc kMu2Sub = &Calculator::getMu2sub;               //
Calculator::ObsFunc kK2NormSub = &Calculator::getK2NormSub;         //
Calculator::ObsFunc kK3Sub = &Calculator::getK3sub;                 //
Calculator::ObsFunc kK3NormSub = &Calculator::getK3NormSub;         //
Calculator::ObsFunc kK3_3sub = &Calculator::getK3_3sub;             //
Calculator::ObsFunc kK3Norm_3sub = &Calculator::getK3Norm_3sub;     //
Calculator::ObsFunc kK3ASub = &Calculator::getK3subA;               //
Calculator::ObsFunc kK3BSub = &Calculator::getK3subB;               //
Calculator::ObsFunc kMu3Sub = &Calculator::getMu3sub;               //
Calculator::ObsFunc kK4Sub = &Calculator::getK4sub;                 //
Calculator::ObsFunc kK4_4sub = &Calculator::getK4_4sub;             //
Calculator::ObsFunc kMu4Sub = &Calculator::getMu4sub;               //
Calculator::ObsFunc kStdSkewSub = &Calculator::getStdSkewSub;       //
Calculator::ObsFunc kIntSkewSub = &Calculator::getIntSkewSub;       //
Calculator::ObsFunc kIntSkewAvgSub = &Calculator::getIntSkewAvgSub; //
Calculator::ObsFunc kExKurSub = &Calculator::getExKurSub;           //
Calculator::ObsFunc kStdSkew_3sub = &Calculator::getStdSkew_3sub;   //
Calculator::ObsFunc kIntSkew_3sub = &Calculator::getIntSkew_3sub;   //
Calculator::ObsFunc kExKur_4sub = &Calculator::getExKur_4sub;       //

// Flow-transverse momentum observables
Calculator::ObsFunc kCV12pt = &Calculator::getCV12pt;   // v1^2-dpt cumulant
Calculator::ObsFunc kCV22pt = &Calculator::getCV22pt;   // v2^2-dpt cumulant
Calculator::ObsFunc kCV32pt = &Calculator::getCV32pt;   // v3^2-dpt cumulant
Calculator::ObsFunc kCV42pt = &Calculator::getCV42pt;   // v4^2-dpt cumulant
Calculator::ObsFunc kCV22pt2 = &Calculator::getCV22pt2; // v2^2-dpt^2 cumulant
Calculator::ObsFunc kCV22pt3 = &Calculator::getCV22pt3; // v2^2-dpt^3 cumulant
Calculator::ObsFunc kCV22pt4 = &Calculator::getCV22pt4; // v2^2-dpt^4 cumulant
Calculator::ObsFunc kCV24pt = &Calculator::getCV24pt;   // v2^4-dpt cumulant
Calculator::ObsFunc kCV24pt2 = &Calculator::getCV24pt2; // v2^4-dpt^2 cumulant

Calculator::ObsFunc kRhoV22pt = &Calculator::getRhoV22pt; // ρ(v2^2,dpt)
Calculator::ObsFunc kRhoV32pt = &Calculator::getRhoV32pt; // ρ(v3^2,dpt)
Calculator::ObsFunc kRhoV42pt = &Calculator::getRhoV42pt; // ρ(v4^2,dpt)

} // namespace GFWSpace
#endif // PWGCF_GENERICFRAMEWORK_CORE_GFWSPACE_H_
