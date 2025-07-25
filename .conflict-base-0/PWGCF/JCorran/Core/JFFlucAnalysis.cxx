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

#include "JFFlucAnalysis.h"

#include <TComplex.h>
#include <TMath.h>

#include <algorithm>

JFFlucAnalysis::JFFlucAnalysis() : TNamed(),
                                   fVertex(0),
                                   fAvgInvariantMass(0.0f),
                                   fCent(0),
                                   fImpactParameter(-1),
                                   subeventMask(kSubEvent_A | kSubEvent_B),
                                   flags(0),
                                   pqvecs(0),
                                   pqvecsRef(0)
{
  //
}

//________________________________________________________________________
JFFlucAnalysis::JFFlucAnalysis(const char* /*name*/) : TNamed(),
                                                       fVertex(0),
                                                       fAvgInvariantMass(0.0f),
                                                       fCent(0),
                                                       fImpactParameter(-1),
                                                       subeventMask(kSubEvent_A | kSubEvent_B),
                                                       flags(0),
                                                       pqvecs(0),
                                                       pqvecsRef(0)
{
  //
}

//________________________________________________________________________
JFFlucAnalysis::JFFlucAnalysis(const JFFlucAnalysis& a) : TNamed(a),
                                                          fVertex(a.fVertex),
                                                          fAvgInvariantMass(a.fAvgInvariantMass),
                                                          fCent(a.fCent),
                                                          fImpactParameter(a.fImpactParameter),
                                                          subeventMask(a.subeventMask),
                                                          flags(a.flags),
                                                          pqvecs(a.pqvecs),
                                                          pqvecsRef(a.pqvecsRef)
{
  // copy constructor
}
//________________________________________________________________________
JFFlucAnalysis& JFFlucAnalysis::operator=(const JFFlucAnalysis& ap)
{
  // assignment operator
  this->~JFFlucAnalysis();
  new (this) JFFlucAnalysis(ap);
  return *this;
}
//________________________________________________________________________
void JFFlucAnalysis::Init()
{
  //
}
//________________________________________________________________________
void JFFlucAnalysis::UserCreateOutputObjects()
{
  //
}

//________________________________________________________________________
JFFlucAnalysis::~JFFlucAnalysis()
{
  //
}

#define C(u) TComplex::Conjugate(u)
// TODO: conjugate macro
inline TComplex TwoGap(const TComplex (&Qa)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], const TComplex (&Qb)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], UInt_t a, UInt_t b)
{
  return Qa[a][1] * C(Qb[b][1]);
}

inline TComplex ThreeGap(const TComplex (&Qa)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], const TComplex (&Qb)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], UInt_t a, UInt_t b, UInt_t c)
{
  return Qa[a][1] * C(Qb[b][1] * Qb[c][1] - Qb[b + c][2]);
}

inline TComplex FourGap22(const TComplex (&Qa)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], const TComplex (&Qb)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], UInt_t a, UInt_t b, UInt_t c, UInt_t d)
{
  return Qa[a][1] * Qa[b][1] * C(Qb[c][1] * Qb[d][1]) - Qa[a + b][2] * C(Qb[c][1] * Qb[d][1]) - Qa[a][1] * Qa[b][1] * C(Qb[c + d][2]) + Qa[a + b][2] * C(Qb[c + d][2]);
}

inline TComplex FourGap13(const TComplex (&Qa)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], const TComplex (&Qb)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], UInt_t a, UInt_t b, UInt_t c, UInt_t d)
{
  return Qa[a][1] * C(Qb[b][1] * Qb[c][1] * Qb[d][1] - Qb[b + c][2] * Qb[d][1] - Qb[b + d][2] * Qb[c][1] - Qb[c + d][2] * Qb[b][1] + 2.0 * Qb[b + c + d][3]);
}

inline TComplex SixGap33(const TComplex (&Qa)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], const TComplex (&Qb)[JFFlucAnalysis::kNH][JFFlucAnalysis::nKL], UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5, UInt_t n6)
{
  return Qa[n1][1] * Qa[n2][1] * Qa[n3][1] * C(Qb[n4][1] * Qb[n5][1] * Qb[n6][1]) - Qa[n1][1] * Qa[n2][1] * Qa[n3][1] * C(Qb[n4 + n5][2] * Qb[n6][1]) - Qa[n1][1] * Qa[n2][1] * Qa[n3][1] * C(Qb[n4 + n6][2] * Qb[n5][1]) - Qa[n1][1] * Qa[n2][1] * Qa[n3][1] * C(Qb[n5 + n6][2] * Qb[n4][1]) + 2.0 * Qa[n1][1] * Qa[n2][1] * Qa[n3][1] * C(Qb[n4 + n5 + n6][3]) - Qa[n1 + n2][2] * Qa[n3][1] * C(Qb[n4][1] * Qb[n5][1] * Qb[n6][1]) + Qa[n1 + n2][2] * Qa[n3][1] * C(Qb[n4 + n5][2] * Qb[n6][1]) + Qa[n1 + n2][2] * Qa[n3][1] * C(Qb[n4 + n6][2] * Qb[n5][1]) + Qa[n1 + n2][2] * Qa[n3][1] * C(Qb[n5 + n6][2] * Qb[n4][1]) - 2.0 * Qa[n1 + n2][2] * Qa[n3][1] * C(Qb[n4 + n5 + n6][3]) - Qa[n1 + n3][2] * Qa[n2][1] * C(Qb[n4][1] * Qb[n5][1] * Qb[n6][1]) + Qa[n1 + n3][2] * Qa[n2][1] * C(Qb[n4 + n5][2] * Qb[n6][1]) + Qa[n1 + n3][2] * Qa[n2][1] * C(Qb[n4 + n6][2] * Qb[n5][1]) + Qa[n1 + n3][2] * Qa[n2][1] * C(Qb[n5 + n6][2] * Qb[n4][1]) - 2.0 * Qa[n1 + n3][2] * Qa[n2][1] * C(Qb[n4 + n5 + n6][3]) - Qa[n2 + n3][2] * Qa[n1][1] * C(Qb[n4][1] * Qb[n5][1] * Qb[n6][1]) + Qa[n2 + n3][2] * Qa[n1][1] * C(Qb[n4 + n5][2] * Qb[n6][1]) + Qa[n2 + n3][2] * Qa[n1][1] * C(Qb[n4 + n6][2] * Qb[n5][1]) + Qa[n2 + n3][2] * Qa[n1][1] * C(Qb[n5 + n6][2] * Qb[n4][1]) - 2.0 * Qa[n2 + n3][2] * Qa[n1][1] * C(Qb[n4 + n5 + n6][3]) + 2.0 * Qa[n1 + n2 + n3][3] * C(Qb[n4][1] * Qb[n5][1] * Qb[n6][1]) - 2.0 * Qa[n1 + n2 + n3][3] * C(Qb[n4 + n5][2] * Qb[n6][1]) - 2.0 * Qa[n1 + n2 + n3][3] * C(Qb[n4 + n6][2] * Qb[n5][1]) - 2.0 * Qa[n1 + n2 + n3][3] * C(Qb[n5 + n6][2] * Qb[n4][1]) + 4.0 * Qa[n1 + n2 + n3][3] * C(Qb[n4 + n5 + n6][3]);
}

TComplex JFFlucAnalysis::Q(int n, int p)
{
  // Return QvectorQC
  // Q{-n, p} = Q{n, p}*
  return n >= 0 ? pqvecs->QvectorQC[n][p] : C(pqvecs->QvectorQC[-n][p]);
}

TComplex JFFlucAnalysis::Q(const JQVectorsT& qvecs, int n, int p)
{
  // Return QvectorQC
  // Q{-n, p} = Q{n, p}*
  return n >= 0 ? qvecs.QvectorQC[n][p] : C(qvecs.QvectorQC[-n][p]);
}

TComplex JFFlucAnalysis::Two(int n1, int n2)
{
  // two-particle correlation <exp[i(n1*phi1 + n2*phi2)]>
  return Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);
}

TComplex JFFlucAnalysis::Four(int n1, int n2, int n3, int n4)
{
  return Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) + Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);
}

TComplex JFFlucAnalysis::TwoDiff(int n1, int n2)
{
#define dp(n, p) Q(*pqvecs, n, p)    // POI
#define dQ(n, p) Q(*pqvecsRef, n, p) // REF
#define dq(n, p) dp(n, p)            //(dp(n,p)+dQ(n,p)) //POI+REF in narrow bin. Since there is no mass for ref, q = POI
                                     // #define dq(n,p) (dp(n,p)+dQ(n,p)) //POI+REF in narrow bin. Since there is no mass for ref, q = POI
  return dp(n1, 1) * dQ(n2, 1) - dq(n1 + n2, 2);
}

TComplex JFFlucAnalysis::FourDiff(int n1, int n2, int n3, int n4)
{
  return dp(n1, 1) * dQ(n2, 1) * dQ(n3, 1) * dQ(n4, 1) - dq(n1 + n2, 2) * dQ(n3, 1) * dQ(n4, 1) - dq(n1 + n3, 2) * dQ(n2, 1) * dQ(n4, 1) - dp(n1, 1) * dQ(n2 + n3, 2) * dQ(n4, 1) + 2. * dq(n1 + n2 + n3, 3) * dQ(n4, 1) - dQ(n2, 1) * dQ(n3, 1) * dq(n1 + n4, 2) + dQ(n2 + n3, 2) * dq(n1 + n4, 2) - dp(n1, 1) * dQ(n3, 1) * dQ(n2 + n4, 2) + dq(n1 + n3, 2) * dQ(n2 + n4, 2) + 2. * dQ(n3, 1) * dq(n1 + n2 + n4, 3) - dp(n1, 1) * dQ(n2, 1) * dQ(n3 + n4, 2) + dq(n1 + n2, 2) * dQ(n3 + n4, 2) + 2. * dQ(n2, 1) * dq(n1 + n3 + n4, 3) + 2. * dp(n1, 1) * dQ(n2 + n3 + n4, 3) - 6. * dq(n1 + n2 + n3 + n4, 4);
}

#undef dp
#undef dQ
#undef dq
#undef C

//________________________________________________________________________
void JFFlucAnalysis::UserExec(Option_t* /*popt*/) // NOLINT(readability/casting) false positive: https://github.com/cpplint/cpplint/issues/131
{
  TComplex corr[kNH][nKL];
  TComplex ncorr[kNH][nKL];
  TComplex ncorr2[kNH][nKL][kcNH][nKL];

  for (UInt_t i = 0; i < 2; ++i) {
    if ((subeventMask & (1 << i)) == 0)
      continue;
    decltype(pqvecs->QvectorQCgap[i])& Qa = pqvecs->QvectorQCgap[i];                                   // this is for one differential bin only.
    decltype(pqvecs->QvectorQCgap[1 - i])& Qb = (pqvecsRef ? pqvecsRef : pqvecs)->QvectorQCgap[1 - i]; // A & B subevents from POI and REF, when given
    Double_t ref_2p = TwoGap(Qa, Qb, 0, 0).Re();
    Double_t ref_3p = ThreeGap(Qa, Qb, 0, 0, 0).Re();
    Double_t ref_4p = FourGap22(Qa, Qb, 0, 0, 0, 0).Re();
    Double_t ref_4pB = FourGap13(Qa, Qb, 0, 0, 0, 0).Re();
    Double_t ref_6p = SixGap33(Qa, Qb, 0, 0, 0, 0, 0, 0).Re();

    Double_t ebe_2p_weight = 1.0;
    Double_t ebe_3p_weight = 1.0;
    Double_t ebe_4p_weight = 1.0;
    Double_t ebe_4p_weightB = 1.0;
    Double_t ebe_6p_weight = 1.0;
    if (flags & kFlucEbEWeighting) {
      ebe_2p_weight = ref_2p;
      ebe_3p_weight = ref_3p;
      ebe_4p_weight = ref_4p;
      ebe_4p_weightB = ref_4pB;
      ebe_6p_weight = ref_6p;
    }
    Double_t ref_2Np[2 * nKL] = {
      ref_2p,
      ref_4p,
      ref_6p};
    Double_t ebe_2Np_weight[2 * nKL] = {
      ebe_2p_weight,
      ebe_4p_weight,
      ebe_6p_weight};
    if (flags & kFlucEbEWeighting) {
      for (UInt_t ik = 3; ik < 2 * nKL; ik++) {
        double dk = static_cast<double>(ik);
        ref_2Np[ik] = ref_2Np[ik - 1] * std::max(Qa[0][1].Re() - dk, 1.0) * std::max(Qb[0][1].Re() - dk, 1.0);
        ebe_2Np_weight[ik] = ebe_2Np_weight[ik - 1] * std::max(Qa[0][1].Re() - dk, 1.0) * std::max(Qb[0][1].Re() - dk, 1.0);
      }
    } else {
      for (UInt_t ik = 3; ik < 2 * nKL; ik++) {
        double dk = static_cast<double>(ik);
        ref_2Np[ik] = ref_2Np[ik - 1] * std::max(Qa[0][1].Re() - dk, 1.0) * std::max(Qb[0][1].Re() - dk, 1.0);
        ebe_2Np_weight[ik] = 1.0;
      }
    }

    for (UInt_t ih = 2; ih < kNH; ih++) {
      corr[ih][1] = TwoGap(Qa, Qb, ih, ih);
      for (UInt_t ik = 2; ik < nKL; ik++)
        corr[ih][ik] = corr[ih][ik - 1] * corr[ih][1]; // TComplex::Power(corr[ih][1],ik);
      ncorr[ih][1] = corr[ih][1];
      ncorr[ih][2] = FourGap22(Qa, Qb, ih, ih, ih, ih);
      ncorr[ih][3] = SixGap33(Qa, Qb, ih, ih, ih, ih, ih, ih);
      for (UInt_t ik = 4; ik < nKL; ik++)
        ncorr[ih][ik] = corr[ih][ik]; // for 8,...-particle correlations, ignore the autocorrelation / weight dependency for now

      for (UInt_t ihh = 2; ihh < kcNH; ihh++) {
        ncorr2[ih][1][ihh][1] = FourGap22(Qa, Qb, ih, ihh, ih, ihh);
        ncorr2[ih][1][ihh][2] = SixGap33(Qa, Qb, ih, ihh, ihh, ih, ihh, ihh);
        ncorr2[ih][2][ihh][1] = SixGap33(Qa, Qb, ih, ih, ihh, ih, ih, ihh);
        for (UInt_t ik = 2; ik < nKL; ik++)
          for (UInt_t ikk = 2; ikk < nKL; ikk++)
            ncorr2[ih][ik][ihh][ikk] = ncorr[ih][ik] * ncorr[ihh][ikk];
      }
    }

    for (UInt_t ih = 2; ih < kNH; ih++) {
      for (UInt_t ik = 1; ik < nKL; ik++) { // 2k(0) =1, 2k(1) =2, 2k(2)=4....
                                            // vn2[ih][ik] = corr[ih][ik].Re() / ref_2Np[ik - 1];
        // fh_vn[ih][ik][fCBin]->Fill(vn2[ih][ik], ebe_2Np_weight[ik - 1]);
        // fh_vna[ih][ik][fCBin]->Fill(ncorr[ih][ik].Re() / ref_2Np[ik - 1], ebe_2Np_weight[ik - 1]);
        phs[HIST_THN_SPARSE_VN]->Fill(fCent, fAvgInvariantMass, ih, ik, ncorr[ih][ik].Re() / ref_2Np[ik - 1], ebe_2Np_weight[ik - 1]);
        for (UInt_t ihh = 2; ihh < kcNH; ihh++) {
          for (UInt_t ikk = 1; ikk < nKL; ikk++) {
            Double_t vn2_vn2 = ncorr2[ih][ik][ihh][ikk] / ref_2Np[ik + ikk - 1];
            phs[HIST_THN_SPARSE_VN_VN]->Fill(fCent, fAvgInvariantMass, ih, ik, ihh, ikk, vn2_vn2, ebe_2Np_weight[ik + ikk - 1]);
          }
        }
      }
    }

    //************************************************************************
    TComplex V4V2star_2 = Qa[4][1] * Qb[2][1] * Qb[2][1];
    TComplex V4V2starv2_2 = V4V2star_2 * corr[2][1] / ref_2Np[0];                           // vn[2][1]
    TComplex V4V2starv2_4 = V4V2star_2 * corr[2][2] / ref_2Np[1];                           // vn2[2][2]
    TComplex V5V2starV3starv2_2 = Qa[5][1] * Qb[2][1] * Qb[3][1] * corr[2][1] / ref_2Np[0]; // vn2[2][1]
    TComplex V5V2starV3star = Qa[5][1] * Qb[2][1] * Qb[3][1];
    TComplex V5V2starV3startv3_2 = V5V2starV3star * corr[3][1] / ref_2Np[0]; // vn2[3][1]
    TComplex V6V2star_3 = Qa[6][1] * Qb[2][1] * Qb[2][1] * Qb[2][1];
    TComplex V6V3star_2 = Qa[6][1] * Qb[3][1] * Qb[3][1];
    TComplex V6V2starV4star = Qa[6][1] * Qb[2][1] * Qb[4][1];
    TComplex V7V2star_2V3star = Qa[7][1] * Qb[2][1] * Qb[2][1] * Qb[3][1];
    TComplex V7V2starV5star = Qa[7][1] * Qb[2][1] * Qb[5][1];
    TComplex V7V3starV4star = Qa[7][1] * Qb[3][1] * Qb[4][1];
    TComplex V8V2starV3star_2 = Qa[8][1] * Qb[2][1] * Qb[3][1] * Qb[3][1];
    TComplex V8V2star_4 = Qa[8][1] * TComplex::Power(Qb[2][1], 4);

    // New correlators (Modified by You's correction term for self-correlations)
    TComplex nV4V2star_2 = ThreeGap(Qa, Qb, 4, 2, 2) / ref_3p;
    TComplex nV5V2starV3star = ThreeGap(Qa, Qb, 5, 2, 3) / ref_3p;
    TComplex nV6V2star_3 = FourGap13(Qa, Qb, 6, 2, 2, 2) / ref_4pB;
    TComplex nV6V3star_2 = ThreeGap(Qa, Qb, 6, 3, 3) / ref_3p;
    TComplex nV6V2starV4star = ThreeGap(Qa, Qb, 6, 2, 4) / ref_3p;
    TComplex nV7V2star_2V3star = FourGap13(Qa, Qb, 7, 2, 2, 3) / ref_4pB;
    TComplex nV7V2starV5star = ThreeGap(Qa, Qb, 7, 2, 5) / ref_3p;
    TComplex nV7V3starV4star = ThreeGap(Qa, Qb, 7, 3, 4) / ref_3p;
    TComplex nV8V2starV3star_2 = FourGap13(Qa, Qb, 8, 2, 3, 3) / ref_4pB;

    TComplex nV4V4V2V2 = FourGap22(Qa, Qb, 4, 2, 4, 2) / ref_4p;
    TComplex nV3V3V2V2 = FourGap22(Qa, Qb, 3, 2, 3, 2) / ref_4p;
    TComplex nV5V5V2V2 = FourGap22(Qa, Qb, 5, 2, 5, 2) / ref_4p;
    TComplex nV5V5V3V3 = FourGap22(Qa, Qb, 5, 3, 5, 3) / ref_4p;
    TComplex nV4V4V3V3 = FourGap22(Qa, Qb, 4, 3, 4, 3) / ref_4p;

    pht[HIST_THN_V4V2starv2_2]->Fill(fCent, fAvgInvariantMass, V4V2starv2_2.Re());
    pht[HIST_THN_V4V2starv2_4]->Fill(fCent, fAvgInvariantMass, V4V2starv2_4.Re());
    pht[HIST_THN_V4V2star_2]->Fill(fCent, fAvgInvariantMass, V4V2star_2.Re(), ebe_3p_weight); // added 2015.3.18
    pht[HIST_THN_V5V2starV3starv2_2]->Fill(fCent, fAvgInvariantMass, V5V2starV3starv2_2.Re());
    pht[HIST_THN_V5V2starV3star]->Fill(fCent, fAvgInvariantMass, V5V2starV3star.Re(), ebe_3p_weight);
    pht[HIST_THN_V5V2starV3startv3_2]->Fill(fCent, fAvgInvariantMass, V5V2starV3startv3_2.Re());
    pht[HIST_THN_V6V2star_3]->Fill(fCent, fAvgInvariantMass, V6V2star_3.Re(), ebe_4p_weightB);
    pht[HIST_THN_V6V3star_2]->Fill(fCent, fAvgInvariantMass, V6V3star_2.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2star_2V3star]->Fill(fCent, fAvgInvariantMass, V7V2star_2V3star.Re(), ebe_4p_weightB);

    pht[HIST_THN_V4V2star_2]->Fill(fCent, fAvgInvariantMass, nV4V2star_2.Re(), ebe_3p_weight); // added 2015.6.10
    pht[HIST_THN_V5V2starV3star]->Fill(fCent, fAvgInvariantMass, nV5V2starV3star.Re(), ebe_3p_weight);
    pht[HIST_THN_V6V3star_2]->Fill(fCent, fAvgInvariantMass, nV6V3star_2.Re(), ebe_3p_weight);

    // use this to avoid self-correlation 4p correlation (2 particles from A, 2 particles from B) -> MA(MA-1)MB(MB-1) : evt weight..
    pht[HIST_THN_nV4V4V2V2]->Fill(fCent, fAvgInvariantMass, nV4V4V2V2.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV3V3V2V2]->Fill(fCent, fAvgInvariantMass, nV3V3V2V2.Re(), ebe_2Np_weight[1]);

    pht[HIST_THN_nV5V5V2V2]->Fill(fCent, fAvgInvariantMass, nV5V5V2V2.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV5V5V3V3]->Fill(fCent, fAvgInvariantMass, nV5V5V3V3.Re(), ebe_2Np_weight[1]);
    pht[HIST_THN_nV4V4V3V3]->Fill(fCent, fAvgInvariantMass, nV4V4V3V3.Re(), ebe_2Np_weight[1]);

    // higher order correlators, added 2017.8.10
    pht[HIST_THN_V8V2starV3star_2]->Fill(fCent, fAvgInvariantMass, V8V2starV3star_2.Re(), ebe_4p_weightB);
    pht[HIST_THN_V8V2star_4]->Fill(fCent, fAvgInvariantMass, V8V2star_4.Re()); // 5p weight
    pht[HIST_THN_V6V2star_3]->Fill(fCent, fAvgInvariantMass, nV6V2star_3.Re(), ebe_4p_weightB);
    pht[HIST_THN_V7V2star_2V3star]->Fill(fCent, fAvgInvariantMass, nV7V2star_2V3star.Re(), ebe_4p_weightB);
    pht[HIST_THN_V8V2starV3star_2]->Fill(fCent, fAvgInvariantMass, nV8V2starV3star_2.Re(), ebe_4p_weightB);

    pht[HIST_THN_V6V2starV4star]->Fill(fCent, fAvgInvariantMass, V6V2starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2starV5star]->Fill(fCent, fAvgInvariantMass, V7V2starV5star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V3starV4star]->Fill(fCent, fAvgInvariantMass, V7V3starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V6V2starV4star]->Fill(fCent, fAvgInvariantMass, nV6V2starV4star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V2starV5star]->Fill(fCent, fAvgInvariantMass, nV7V2starV5star.Re(), ebe_3p_weight);
    pht[HIST_THN_V7V3starV4star]->Fill(fCent, fAvgInvariantMass, nV7V3starV4star.Re(), ebe_3p_weight);

    Double_t event_weight_two_gap = 1.0;
    if (flags & kFlucEbEWeighting) {
      event_weight_two_gap = (Qa[0][1] * Qb[0][1]).Re();
    }

    for (UInt_t ih = 2; ih < kNH; ih++) {
      TComplex sctwoGap = (Qa[ih][1] * TComplex::Conjugate(Qb[ih][1])) / (Qa[0][1] * Qb[0][1]).Re();
      pht[HIST_THN_SC_with_QC_2corr_gap]->Fill(fCent, fAvgInvariantMass, ih, sctwoGap.Re(), event_weight_two_gap);
    }
  }

  auto four = [&](int a, int b, int c, int d) -> TComplex { return pqvecsRef ? FourDiff(a, b, c, d) : Four(a, b, c, d); };
  auto two = [&](int a, int b) -> TComplex { return pqvecsRef ? TwoDiff(a, b) : Two(a, b); };
  Double_t event_weight_four = 1.0;
  Double_t event_weight_two = 1.0;
  if (flags & kFlucEbEWeighting) {
    event_weight_four = four(0, 0, 0, 0).Re();
    event_weight_two = two(0, 0).Re();
  }

  for (UInt_t ih = 2; ih < kNH; ih++) {
    for (UInt_t ihh = 2, mm = (ih < kcNH ? ih : static_cast<UInt_t>(kcNH)); ihh < mm; ihh++) {
      TComplex scfour = four(ih, ihh, -ih, -ihh) / four(0, 0, 0, 0).Re();
      pht[HIST_THN_SC_with_QC_4corr]->Fill(fCent, fAvgInvariantMass, ih, ihh, scfour.Re(), event_weight_four);
    }
    TComplex sctwo = two(ih, -ih) / two(0, 0).Re();
    pht[HIST_THN_SC_with_QC_2corr]->Fill(fCent, fAvgInvariantMass, ih, sctwo.Re(), event_weight_two);
  }
}

//________________________________________________________________________
void JFFlucAnalysis::Terminate(Option_t* /*popt*/) // NOLINT(readability/casting) false positive: https://github.com/cpplint/cpplint/issues/131
{
  //
}
