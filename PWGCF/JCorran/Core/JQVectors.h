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

#ifndef PWGCF_JCORRAN_CORE_JQVECTORS_H_
#define PWGCF_JCORRAN_CORE_JQVECTORS_H_

#include <experimental/type_traits>
#include <TMath.h>

template <class Q, UInt_t nh, UInt_t nk>
class JQVectorsGapBase
{
 public:
  JQVectorsGapBase() {}
  virtual ~JQVectorsGapBase() {}
  Q QvectorQCgap[2][nh][nk];
};

class JQVectorsEmptyBase
{
 public:
  //
};

template <class Q, UInt_t nh, UInt_t nk, bool gap>
class JQVectors : public std::conditional_t<gap, JQVectorsGapBase<Q, nh, nk>, JQVectorsEmptyBase>
{
 public:
  JQVectors() {}
  ~JQVectors() {}

  template <class T>
  using hasWeightNUA = decltype(std::declval<T&>().weightNUA());
  template <class T>
  using hasWeightEff = decltype(std::declval<T&>().weightEff());
  template <class T>
  using hasInvMass = decltype(std::declval<T&>().invMass());

  template <class JInputClass>
  inline void Calculate(JInputClass& inputInst, float etamin, float etamax, float massMin = 0.0f, float massMax = 999.9f)
  {
    // calculate Q-vector for QC method ( no subgroup )
    for (UInt_t ih = 0; ih < nh; ++ih) {
      for (UInt_t ik = 0; ik < nk; ++ik) {
        QvectorQC[ih][ik] = Q(0, 0);
        if constexpr (gap) {
          for (UInt_t isub = 0; isub < 2; ++isub)
            this->QvectorQCgap[isub][ih][ik] = Q(0, 0);
        }
      }
    }
    for (auto& track : inputInst) {
      if (track.eta() < -etamax || track.eta() > etamax)
        continue;
      using JInputClassIter = typename JInputClass::iterator;
      if constexpr (std::experimental::is_detected<hasInvMass, const JInputClassIter>::value) {
        if (track.invMass() < massMin || track.invMass() >= massMax)
          continue;
      }

      UInt_t isub = (UInt_t)(track.eta() > 0.0);
      for (UInt_t ih = 0; ih < nh; ++ih) {
        Double_t tf = 1.0;
        for (UInt_t ik = 0; ik < nk; ++ik) {
          Q q(tf * TMath::Cos(ih * track.phi()), tf * TMath::Sin(ih * track.phi()));
          QvectorQC[ih][ik] += q;

          if constexpr (gap) {
            if (TMath::Abs(track.eta()) > etamin)
              this->QvectorQCgap[isub][ih][ik] += q;
          }

          if constexpr (std::experimental::is_detected<hasWeightNUA, const JInputClassIter>::value)
            tf /= track.weightNUA();
          if constexpr (std::experimental::is_detected<hasWeightEff, const JInputClassIter>::value)
            tf *= track.weightEff();
        }
      }
    }
  }
  Q QvectorQC[nh][nk];
};

#endif // PWGCF_JCORRAN_CORE_JQVECTORS_H_
