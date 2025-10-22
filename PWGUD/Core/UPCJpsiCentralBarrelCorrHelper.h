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
///
/// \brief
/// \author Sara Haidlova, sara.haidlova@cern.ch
/// \since March 2024

#ifndef PWGUD_CORE_UPCJPSICENTRALBARRELCORRHELPER_H_
#define PWGUD_CORE_UPCJPSICENTRALBARRELCORRHELPER_H_

#include "CommonConstants/MathConstants.h"

#include <algorithm>
#include <random>
#include <vector>

/*enum ParticleType {
  P_ELECTRON = 0,
  P_MUON = 1,
  P_PROTON = 2
};*/

/*template <typename T>
int testPIDhypoTPC(T trackPID)
{
  float nSigmaTPC[3];
  nSigmaTPC[P_ELECTRON] = std::abs(trackPID.tpcNSigmaEl());
  nSigmaTPC[P_MUON] = std::abs(trackPID.tpcNSigmaMu());
  nSigmaTPC[P_PROTON] = std::abs(trackPID.tpcNSigmaPr());
  int enumChoiceTPC = std::distance(std::begin(nSigmaTPC),
                                    std::min_element(std::begin(nSigmaTPC), std::end(nSigmaTPC)));
  if (trackPID.hasTPC()) {
    return enumChoiceTPC;
  } else {
    return -1;
  }
}

template <typename T>
int testPIDhypo(T trackPID)
{
  float nSigmaTPC[3];
  nSigmaTPC[P_ELECTRON] = std::abs(trackPID.tpcNSigmaEl());
  nSigmaTPC[P_MUON] = std::abs(trackPID.tpcNSigmaMu());
  nSigmaTPC[P_PROTON] = std::abs(trackPID.tpcNSigmaPr());
  int enumChoiceTPC = std::distance(std::begin(nSigmaTPC),
                                    std::min_element(std::begin(nSigmaTPC), std::end(nSigmaTPC)));

  float nSigmaTOF[3];
  nSigmaTOF[P_ELECTRON] = std::abs(trackPID.tofNSigmaEl());
  nSigmaTOF[P_MUON] = std::abs(trackPID.tofNSigmaMu());
  nSigmaTOF[P_PROTON] = std::abs(trackPID.tofNSigmaPr());
  int enumChoiceTOF = std::distance(std::begin(nSigmaTOF),
                                    std::min_element(std::begin(nSigmaTOF), std::end(nSigmaTOF)));
  if (trackPID.hasTPC() || trackPID.hasTOF()) {
    if (trackPID.hasTOF()) {
      return enumChoiceTOF;
    } else {
      return enumChoiceTPC;
    }
  } else {
    return -1;
  }
}*/

float* correlation(TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
{
  TLorentzVector pa(1., 0., 0, 1);  // projectile
  TLorentzVector pb(1., 0., 0, -1); // target

  float* q = new float[3]; // to save values

  // Accoplanarity angle
  Float_t deltaPhi;
  deltaPhi = lv1->Phi() - lv2->Phi();
  float accoplCut = 1. - std::abs(deltaPhi) / o2::constants::math::PI;
  // z
  TLorentzVector z;
  Float_t part1, part2;

  // Dot product: v1*v2 = t1*t2-x1*x2-y1*y2-z1*z2

  part1 = lv->Dot(pb);
  part2 = lv->Dot(pa);

  Float_t part3x = pa.X() * part1;
  Float_t part3y = pa.Y() * part1;
  Float_t part3z = pa.Z() * part1;
  Float_t part3e = pa.T() * part1;

  Float_t part4x = pb.X() * part2;
  Float_t part4y = pb.Y() * part2;
  Float_t part4z = pb.Z() * part2;
  Float_t part4e = pb.T() * part2;

  TLorentzVector part3(TVector3(part3x, part3y, part3z), part3e);
  TLorentzVector part4(TVector3(part4x, part4y, part4z), part4e);

  // Un-normalized Z
  z = part3 - part4;

  // Normalized z
  Float_t normz = std::sqrt(-z * z);
  Float_t znx = z.X() / normz;
  Float_t zny = z.Y() / normz;
  Float_t znz = z.Z() / normz;
  Float_t zne = z.E() / normz;

  // Normalized z
  TLorentzVector zhat(TVector3(znx, zny, znz), zne);

  // calculate x
  TLorentzVector x;

  Float_t constant1 = (lv->Dot(*lv)) / (2 * (lv->Dot(pa)));
  Float_t constant2 = (lv->Dot(*lv)) / (2 * (lv->Dot(pb)));

  Float_t comp1x = pa.X() * constant1;
  Float_t comp1y = pa.Y() * constant1;
  Float_t comp1z = pa.Z() * constant1;
  Float_t comp1e = pa.T() * constant1;

  TLorentzVector comp1(TVector3(comp1x, comp1y, comp1z), comp1e);

  Float_t comp2x = pb.X() * constant2;
  Float_t comp2y = pb.Y() * constant2;
  Float_t comp2z = pb.Z() * constant2;
  Float_t comp2e = pb.T() * constant2;

  TLorentzVector comp2(TVector3(comp2x, comp2y, comp2z), comp2e);

  // Un-normalized x
  x = *lv - comp1 - comp2;
  // normalize x
  Float_t normx = std::sqrt(-x * x);
  Float_t xnx = x.X() / normx;
  Float_t xny = x.Y() / normx;
  Float_t xnz = x.Z() / normx;
  Float_t xne = x.E() / normx;

  // Normalized x
  TLorentzVector xhat(TVector3(xnx, xny, xnz), xne);

  // calculate y
  // TLorentzVector y;
  Float_t yone = pa.Y() * pb.Z() * lv->E() - pa.Z() * pb.Y() * lv->E() + pa.Z() * pb.E() * lv->Y() + pa.E() * pb.Y() * lv->Z() - pa.Y() * pb.E() * lv->Z() - pa.E() * pb.Z() * lv->Y();
  Float_t ytwo = -pa.Z() * pb.E() * lv->X() + pa.Z() * pb.X() * lv->E() - pa.X() * pb.Z() * lv->E() + pa.X() * pb.E() * lv->Z() - pa.E() * pb.X() * lv->Z() + pa.E() * pb.Z() * lv->X();
  Float_t ythree = pa.X() * pb.Y() * lv->E() - pa.Y() * pb.X() * lv->E() + pa.Y() * pb.E() * lv->X() - pa.X() * pb.E() * lv->Y() + pa.E() * pb.X() * lv->Y() - pa.E() * pb.Y() * lv->X();
  Float_t yfour = -pa.X() * pb.Y() * lv->Z() + pa.X() * pb.Z() * lv->Y() - pa.Z() * pb.X() * lv->Y() + pa.Z() * pb.Y() * lv->X() - pa.Y() * pb.Z() * lv->X() + pa.Y() * pb.X() * lv->Z();

  // Un-normalized y
  TLorentzVector y(TVector3(yone, ytwo, ythree), yfour);

  // normalize y
  Float_t normy = std::sqrt(-y * y);
  Float_t ynx = y.X() / normy;
  Float_t yny = y.Y() / normy;
  Float_t ynz = y.Z() / normy;
  Float_t yne = y.E() / normy;

  // normalized y
  TLorentzVector yhat(TVector3(ynx, yny, ynz), yne);

  // Lepton momentum difference
  TLorentzVector diff;
  diff = (*lv1 - *lv2);
  Float_t diff2x = diff.X() / 2.;
  Float_t diff2y = diff.Y() / 2.;
  Float_t diff2z = diff.Z() / 2.;
  Float_t diff2e = diff.E() / 2.;
  TLorentzVector diff2(TVector3(diff2x, diff2y, diff2z), diff2e);

  // Normalize diff2
  Float_t norm2 = std::sqrt(-diff2 * diff2);
  Float_t diff3x = diff2.X() / norm2;
  Float_t diff3y = diff2.Y() / norm2;
  Float_t diff3z = diff2.Z() / norm2;
  Float_t diff3e = diff2.E() / norm2;

  TLorentzVector diff3(TVector3(diff3x, diff3y, diff3z), diff3e);

  // computing the angles
  float cosThetaCS = zhat * diff3;
  Double_t SinThetaCosPhiCS = xhat * diff3;
  Double_t SinThetaSinPhiCS = yhat * diff3;
  //**************************************

  float phi = atan2(SinThetaSinPhiCS, SinThetaCosPhiCS);
  // if (phi>=0) phi = phi;
  if (phi < 0)
    phi = phi + o2::constants::math::TwoPI;

  q[0] = accoplCut;
  q[1] = phi;
  q[2] = cosThetaCS;

  return q;
}

double DeltaPhi(TLorentzVector lv1, TLorentzVector lv2)
{
  TLorentzVector lv_sum = lv1 + lv2;
  TLorentzVector lv_diff = lv1 - lv2;

  double dp = lv_sum.DeltaPhi(lv_diff);

  return dp;
}

double DeltaPhiRandom(TLorentzVector lv1, TLorentzVector lv2)
{
  std::vector<int> indices = {0, 1};
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));
  std::array<TLorentzVector, 2> arrayLorentz = {lv1, lv2};
  TLorentzVector lv_sum = arrayLorentz[indices[0]] + arrayLorentz[indices[1]];
  TLorentzVector lv_diff = arrayLorentz[indices[0]] - arrayLorentz[indices[1]];
  ;

  double dp = lv_sum.DeltaPhi(lv_diff);

  return dp;
}

#endif // PWGUD_CORE_UPCJPSICENTRALBARRELCORRHELPER_H_
