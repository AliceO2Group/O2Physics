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
//
/// \brief utility functions for femto task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 30 May 2023

#ifndef PWGCF_FEMTO3D_CORE_FEMTO3DPAIRTASK_H_
#define PWGCF_FEMTO3D_CORE_FEMTO3DPAIRTASK_H_

// #include "Framework/ASoA.h"
// #include "Framework/DataTypes.h"
// #include "Framework/AnalysisDataModel.h"
// #include "Common/DataModel/PIDResponse.h"
// #include "Framework/Logger.h"
// #include "Common/DataModel/Multiplicity.h"

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"

double particle_mass(int PDGcode)
{
  // if(PDGcode == 2212) return TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (PDGcode == 1000010020)
    return 1.87561294257;
  else
    return TDatabasePDG::Instance()->GetParticle(PDGcode)->Mass();
}

namespace o2::aod::singletrackselector
{
template <typename Type>
Type getBinIndex(float const& value, std::vector<float> const& binning, int const& NsubBins = 1)
{
  Type res = -100;
  for (int i = 0; i < binning.size() - 1; i++) {
    if (value >= binning[i] && binning[i + 1] > value) {
      if (NsubBins < 2) {
        res = (Type)i;
        break;
      } else {
        float subBinWidth = (binning[i + 1] - binning[i]) / NsubBins;
        int subBin = std::floor((value - binning[i]) / subBinWidth);
        int delimeter = std::pow(10, std::to_string(NsubBins).size());

        res = (Type)i + (Type)subBin / delimeter;
        break;
      }
    }
  }
  return res;
}

//====================================================================================

float GetKstarFrom4vectors(TLorentzVector& first4momentum, TLorentzVector& second4momentum, bool isIdentical)
{
  if (isIdentical) {
    TLorentzVector fourmomentadiff = first4momentum - second4momentum;
    return 0.5 * abs(fourmomentadiff.Mag());
  } else {
    TLorentzVector fourmomentasum = first4momentum + second4momentum;
    TLorentzVector fourmomentadif = first4momentum - second4momentum;

    fourmomentadif.Boost((-1) * fourmomentasum.BoostVector());

    return 0.5 * abs(fourmomentadif.Vect().Mag());
  }
}

//====================================================================================

TVector3 GetQLCMSFrom4vectors(TLorentzVector& first4momentum, TLorentzVector& second4momentum)
{
  TLorentzVector fourmomentasum = first4momentum + second4momentum;
  TLorentzVector fourmomentadif = first4momentum - second4momentum;

  fourmomentadif.Boost(0.0, 0.0, (-1) * fourmomentasum.BoostVector().Z()); // boost to LCMS
  fourmomentadif.RotateZ((-1) * fourmomentasum.Phi());                     // rotate so the X axis is along pair's kT

  return fourmomentadif.Vect();
}

//====================================================================================

template <typename TrackType>
class FemtoPair
{
 public:
  FemtoPair() {}
  FemtoPair(TrackType const& first, TrackType const& second)
  {
    _first = first;
    _second = second;
  }
  FemtoPair(TrackType const& first, TrackType const& second, const bool& isidentical)
  {
    _first = first;
    _second = second;
    _isidentical = isidentical;
  }

  FemtoPair(const FemtoPair& obj)
  {
    SetFirstParticle(obj.GetFirstParticle());
    SetSecondParticle(obj.GetSecondParticle());
  }
  explicit FemtoPair(const FemtoPair* obj)
  {
    SetFirstParticle(obj->GetFirstParticle());
    SetSecondParticle(obj->GetSecondParticle());
  }
  ~FemtoPair() {}
  FemtoPair& operator=(const FemtoPair& obj)
  {
    if (this != &obj) {
      SetFirstParticle(obj.GetFirstParticle());
      SetSecondParticle(obj.GetSecondParticle());
    }
    return *this;
  }

  void SetPair(TrackType const& first, TrackType const& second)
  {
    _first = first;
    _second = second;
  }
  void SetFirstParticle(TrackType const& first) { _first = first; }
  void SetSecondParticle(TrackType const& second) { _second = second; }
  void SetIdentical(const bool& isidentical) { _isidentical = isidentical; }
  void SetMagField1(const float& magfield1) { _magfield1 = magfield1; }
  void SetMagField2(const float& magfield2) { _magfield2 = magfield2; }
  void SetPDG1(const int& PDG1) { _PDG1 = PDG1; }
  void SetPDG2(const int& PDG2) { _PDG2 = PDG2; }
  int GetPDG1() { return _PDG1; }
  int GetPDG2() { return _PDG2; }
  void ResetPair();
  void ResetAll();

  TrackType* GetFirstParticle() const { return _first; }
  TrackType* GetSecondParticle() const { return _second; }
  bool IsIdentical() { return _isidentical; }

  bool IsClosePair(const float& deta = 0.01, const float& dphi = 0.01, const float& radius = 1.2) const;
  float GetEtaDiff() const
  {
    if (_first != NULL && _second != NULL)
      return _first->eta() - _second->eta();
    else
      return 1000;
  }
  float GetPhiStarDiff(const float& radius = 1.2) const
  {
    if (_first != NULL && _second != NULL)
      return _first->phiStar(_magfield1, radius) - _second->phiStar(_magfield2, radius);
    else
      return 1000;
  }
  float GetKstar() const;
  TVector3 GetQLCMS() const;
  float GetKt() const;
  float GetMt() const; // test

 private:
  TrackType _first = NULL;
  TrackType _second = NULL;
  float _magfield1 = 0.0, _magfield2 = 0.0;
  int _PDG1 = 0, _PDG2 = 0;
  bool _isidentical = true;
};

template <typename TrackType>
void FemtoPair<TrackType>::ResetPair()
{
  _first = NULL;
  _second = NULL;
}

template <typename TrackType>
void FemtoPair<TrackType>::ResetAll()
{
  _first = NULL;
  _second = NULL;
  _magfield1 = 0.0;
  _magfield2 = 0.0;
  _PDG1 = 0;
  _PDG2 = 0;
  _isidentical = true;
}

template <typename TrackType>
bool FemtoPair<TrackType>::IsClosePair(const float& deta, const float& dphi, const float& radius) const
{
  if (_first == NULL || _second == NULL)
    return true;
  if (_magfield1 * _magfield2 == 0)
    return true;
  if (abs(GetEtaDiff()) < deta && abs(GetPhiStarDiff(radius)) < dphi)
    return true;

  return false;
}

template <typename TrackType>
float FemtoPair<TrackType>::GetKstar() const
{
  if (_first == NULL || _second == NULL)
    return -1000;
  if (_magfield1 * _magfield2 == 0)
    return -1000;
  if (_PDG1 * _PDG2 == 0)
    return -1000;

  TLorentzVector first4momentum;
  first4momentum.SetPtEtaPhiM(_first->pt(), _first->eta(), _first->phi(), particle_mass(_PDG1));
  TLorentzVector second4momentum;
  second4momentum.SetPtEtaPhiM(_second->pt(), _second->eta(), _second->phi(), particle_mass(_PDG2));

  return GetKstarFrom4vectors(first4momentum, second4momentum, _isidentical);
}

template <typename TrackType>
TVector3 FemtoPair<TrackType>::GetQLCMS() const
{
  if (_first == NULL || _second == NULL)
    return TVector3(-1000, -1000, -1000);
  if (_magfield1 * _magfield2 == 0)
    return TVector3(-1000, -1000, -1000);
  if (_PDG1 * _PDG2 == 0)
    return TVector3(-1000, -1000, -1000);

  TLorentzVector first4momentum;
  first4momentum.SetPtEtaPhiM(_first->pt(), _first->eta(), _first->phi(), particle_mass(_PDG1));
  TLorentzVector second4momentum;
  second4momentum.SetPtEtaPhiM(_second->pt(), _second->eta(), _second->phi(), particle_mass(_PDG2));

  return GetQLCMSFrom4vectors(first4momentum, second4momentum);
}

template <typename TrackType>
float FemtoPair<TrackType>::GetKt() const
{
  if (_first == NULL || _second == NULL)
    return -1000;
  if (_magfield1 * _magfield2 == 0)
    return -1000;
  if (_PDG1 * _PDG2 == 0)
    return -1000;

  return 0.5 * std::sqrt(std::pow(_first->px() + _second->px(), 2) + std::pow(_first->py() + _second->py(), 2));
}

template <typename TrackType>
float FemtoPair<TrackType>::GetMt() const
{
  if (_first == NULL || _second == NULL)
    return -1000;
  if (_magfield1 * _magfield2 == 0)
    return -1000;
  if (_PDG1 * _PDG2 == 0)
    return -1000;

  TLorentzVector first4momentum;
  first4momentum.SetPtEtaPhiM(_first->pt(), _first->eta(), _first->phi(), particle_mass(_PDG1));
  TLorentzVector second4momentum;
  second4momentum.SetPtEtaPhiM(_second->pt(), _second->eta(), _second->phi(), particle_mass(_PDG2));

  TLorentzVector fourmomentasum = first4momentum + second4momentum;

  return 0.5 * fourmomentasum.Mt();
}
} // namespace o2::aod::singletrackselector

#endif // PWGCF_FEMTO3D_CORE_FEMTO3DPAIRTASK_H_
