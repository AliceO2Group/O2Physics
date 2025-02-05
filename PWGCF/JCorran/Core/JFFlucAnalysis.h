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
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \since Sep 2022

#ifndef PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_
#define PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_

#include <experimental/type_traits>
#include "JQVectors.h"
#include <TComplex.h>
#include <TNamed.h>
#include <TH1.h>
#include <THn.h>
#include <THnSparse.h>

class JFFlucAnalysis : public TNamed
{
 protected:
  JFFlucAnalysis();
  explicit JFFlucAnalysis(const char* name);
  explicit JFFlucAnalysis(const JFFlucAnalysis& a);    // not implemented
  JFFlucAnalysis& operator=(const JFFlucAnalysis& ap); // not implemented
 public:
  ~JFFlucAnalysis();
  void UserCreateOutputObjects();
  void Init();
  TComplex Q(int n, int p);
  TComplex Two(int n1, int n2);
  TComplex Four(int n1, int n2, int n3, int n4);
  TComplex TwoDiff(int n1, int n2);
  TComplex FourDiff(int n1, int n2, int n3, int n4);
  void UserExec(Option_t* option);
  void Terminate(Option_t*);

  inline void SetEventCentrality(float cent) { fCent = cent; }
  inline float GetEventCentrality() const { return fCent; }
  inline void SetEventImpactParameter(float ip) { fImpactParameter = ip; }
  inline void SetEventVertex(float zvertex) { fVertex = zvertex; }
  inline void SetAverageInvariantMass(float mass) { fAvgInvariantMass = mass; }
  enum SubEvent {
    kSubEvent_A = 0x1,
    kSubEvent_B = 0x2
  };
  inline void SelectSubevents(UInt_t _subeventMask) { subeventMask = _subeventMask; }
  enum HIST_TH1 {
    HIST_TH1_CENTRALITY,
    HIST_TH1_IMPACTPARAM,
    HIST_TH1_ZVERTEX,
    HIST_TH1_COUNT
  };
  enum HIST_THN {
    HIST_THN_PHIETAZ,
    HIST_THN_PTETA,
    HIST_THN_PHIETA,
    HIST_THN_SC_with_QC_4corr,
    HIST_THN_SC_with_QC_2corr,
    HIST_THN_SC_with_QC_2corr_gap,
    HIST_THN_V4V2star_2,
    HIST_THN_V4V2starv2_2,
    HIST_THN_V4V2starv2_4,
    HIST_THN_V5V2starV3starv2_2,
    HIST_THN_V5V2starV3star,
    HIST_THN_V5V2starV3startv3_2,
    HIST_THN_V6V2star_3,
    HIST_THN_V6V3star_2,
    HIST_THN_V6V2starV4star,
    HIST_THN_V7V2star_2V3star,
    HIST_THN_V7V2starV5star,
    HIST_THN_V7V3starV4star,
    HIST_THN_V8V2starV3star_2,
    HIST_THN_V8V2star_4,
    HIST_THN_nV4V2star_2,
    HIST_THN_nV5V2starV3star,
    HIST_THN_nV6V2star_3,
    HIST_THN_nV6V3star_2,
    HIST_THN_nV6V2starV4star,
    HIST_THN_nV7V2star_2V3star,
    HIST_THN_nV7V2starV5star,
    HIST_THN_nV7V3starV4star,
    HIST_THN_nV8V2starV3star_2,
    HIST_THN_nV4V4V2V2,
    HIST_THN_nV3V3V2V2,
    HIST_THN_nV5V5V2V2,
    HIST_THN_nV5V5V3V3,
    HIST_THN_nV4V4V3V3,
    HIST_THN_COUNT
  };
  enum HIST_THN_SPARSE {
    HIST_THN_SPARSE_VN,
    HIST_THN_SPARSE_VN_VN,
    HIST_THN_SPARSE_COUNT
  };
  enum {
    kFlucEbEWeighting = 0x1
  };
  inline void AddFlags(UInt_t _flags) { flags |= _flags; }

  enum { kH0,
         kH1,
         kH2,
         kH3,
         kH4,
         kH5,
         kH6,
         kH7,
         kH8,
         kH9,
         kH10,
         kH11,
         kH12,
         kNH }; // harmonics
  enum { kK0,
         kK1,
         kK2,
         kK3,
         kK4,
         nKL }; // order
  using JQVectorsT = JQVectors<TComplex, kNH, nKL, true>;
  TComplex Q(const JQVectorsT& qvecs, int n, int p);
  inline void SetJQVectors(const JQVectorsT* _pqvecs)
  {
    pqvecs = _pqvecs;
    pqvecsRef = 0;
  }
  inline void SetJQVectors(const JQVectorsT* _pqvecs, const JQVectorsT* _pqvecsRef)
  {
    pqvecs = _pqvecs;
    pqvecsRef = _pqvecsRef;
  }

  template <class T>
  using hasWeightNUA = decltype(std::declval<T&>().weightNUA());
  template <class T>
  using hasWeightEff = decltype(std::declval<T&>().weightEff());
  template <class T>
  using hasSign = decltype(std::declval<T&>().sign());

  template <class JInputClass>
  inline void FillQA(JInputClass& inputInst, UInt_t type = 0u)
  {
    ph1[HIST_TH1_CENTRALITY]->Fill(fCent);
    ph1[HIST_TH1_IMPACTPARAM]->Fill(fImpactParameter);

    for (auto& track : inputInst) {
      Double_t corrInv = 1.0;
      using JInputClassIter = typename JInputClass::iterator;
      if constexpr (std::experimental::is_detected<hasWeightEff, const JInputClassIter>::value)
        corrInv /= track.weightEff();
      if constexpr (std::experimental::is_detected<hasSign, const JInputClassIter>::value)
        pht[HIST_THN_PTETA]->Fill(fCent, track.pt(), track.eta(), track.sign(), corrInv);
      else
        pht[HIST_THN_PTETA]->Fill(fCent, track.pt(), track.eta(), 0.0, corrInv);
      if constexpr (std::experimental::is_detected<hasWeightNUA, const JInputClassIter>::value)
        corrInv /= track.weightNUA();
      pht[HIST_THN_PHIETA]->Fill(fCent, track.phi(), track.eta(), corrInv);
      pht[HIST_THN_PHIETAZ]->Fill(fCent, static_cast<Double_t>(type), track.phi(), track.eta(), fVertex, corrInv);
    }

    ph1[HIST_TH1_ZVERTEX]->Fill(fVertex);
  }

#define kcNH kH4 // max second dimension + 1
 protected:
  Float_t fVertex;           //!
  Float_t fAvgInvariantMass; //!
  Float_t fCent;             //!
  Float_t fImpactParameter;  //!
  UInt_t subeventMask;       //!
  UInt_t flags;              //!

  const JQVectorsT* pqvecs;    //!
  const JQVectorsT* pqvecsRef; //!

  TH1* ph1[HIST_TH1_COUNT];              //!
  THnSparse* pht[HIST_THN_COUNT];        //!
  THnSparse* phs[HIST_THN_SPARSE_COUNT]; //!

  ClassDef(JFFlucAnalysis, 1)
};

#endif // PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_
