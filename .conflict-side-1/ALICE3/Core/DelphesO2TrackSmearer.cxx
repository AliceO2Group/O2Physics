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
/// @file DelphesO2TrackSmearer.cxx
/// @brief Porting to O2Physics of DelphesO2 code.
///        Minimal changes have been made to the original code for adaptation purposes, formatting and commented parts have been considered.
///        Relevant sources:
///                 DelphesO2/src/lutCovm.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutCovm.hh
///                 DelphesO2/src/TrackSmearer.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.cc
///                 DelphesO2/src/TrackSmearer.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.hh
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it
///

//////////////////////////////
// DelphesO2/src/lutCovm.cc //
//////////////////////////////

/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

// #include "TrackSmearer.hh"
// #include "TrackUtils.hh"
// #include "TRandom.h"
// #include <iostream>
// #include <fstream>

#include "ALICE3/Core/DelphesO2TrackSmearer.h"

namespace o2
{
namespace delphes
{

/*****************************************************************/

bool TrackSmearer::loadTable(int pdg, const char* filename, bool forceReload)
{
  auto ipdg = getIndexPDG(pdg);
  if (mLUTHeader[ipdg] && !forceReload) {
    std::cout << " --- LUT table for PDG " << pdg << " has been already loaded with index " << ipdg << std::endl;
    return false;
  }
  mLUTHeader[ipdg] = new lutHeader_t;

  std::ifstream lutFile(filename, std::ifstream::binary);
  if (!lutFile.is_open()) {
    std::cout << " --- cannot open covariance matrix file for PDG " << pdg << ": " << filename << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  lutFile.read(reinterpret_cast<char*>(mLUTHeader[ipdg]), sizeof(lutHeader_t));
  if (lutFile.gcount() != sizeof(lutHeader_t)) {
    std::cout << " --- troubles reading covariance matrix header for PDG " << pdg << ": " << filename << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  if (mLUTHeader[ipdg]->version != LUTCOVM_VERSION) {
    std::cout << " --- LUT header version mismatch: expected/detected = " << LUTCOVM_VERSION << "/" << mLUTHeader[ipdg]->version << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  if (mLUTHeader[ipdg]->pdg != pdg) {
    std::cout << " --- LUT header PDG mismatch: expected/detected = " << pdg << "/" << mLUTHeader[ipdg]->pdg << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  const int nnch = mLUTHeader[ipdg]->nchmap.nbins;
  const int nrad = mLUTHeader[ipdg]->radmap.nbins;
  const int neta = mLUTHeader[ipdg]->etamap.nbins;
  const int npt = mLUTHeader[ipdg]->ptmap.nbins;
  mLUTEntry[ipdg] = new lutEntry_t****[nnch];
  for (int inch = 0; inch < nnch; ++inch) {
    mLUTEntry[ipdg][inch] = new lutEntry_t***[nrad];
    for (int irad = 0; irad < nrad; ++irad) {
      mLUTEntry[ipdg][inch][irad] = new lutEntry_t**[neta];
      for (int ieta = 0; ieta < neta; ++ieta) {
        mLUTEntry[ipdg][inch][irad][ieta] = new lutEntry_t*[npt];
        for (int ipt = 0; ipt < npt; ++ipt) {
          mLUTEntry[ipdg][inch][irad][ieta][ipt] = new lutEntry_t;
          lutFile.read(reinterpret_cast<char*>(mLUTEntry[ipdg][inch][irad][ieta][ipt]), sizeof(lutEntry_t));
          if (lutFile.gcount() != sizeof(lutEntry_t)) {
            std::cout << " --- troubles reading covariance matrix entry for PDG " << pdg << ": " << filename << std::endl;
            return false;
          }
        }
      }
    }
  }
  std::cout << " --- read covariance matrix table for PDG " << pdg << ": " << filename << std::endl;
  mLUTHeader[ipdg]->print();

  lutFile.close();
  return true;
}

/*****************************************************************/

lutEntry_t*
  TrackSmearer::getLUTEntry(int pdg, float nch, float radius, float eta, float pt, float& interpolatedEff)
{
  auto ipdg = getIndexPDG(pdg);
  if (!mLUTHeader[ipdg])
    return nullptr;
  auto inch = mLUTHeader[ipdg]->nchmap.find(nch);
  auto irad = mLUTHeader[ipdg]->radmap.find(radius);
  auto ieta = mLUTHeader[ipdg]->etamap.find(eta);
  auto ipt = mLUTHeader[ipdg]->ptmap.find(pt);

  // Interpolate if requested
  auto fraction = mLUTHeader[ipdg]->nchmap.fracPositionWithinBin(nch);
  if (mInterpolateEfficiency) {
    if (fraction > 0.5) {
      if (mWhatEfficiency == 1) {
        if (inch < mLUTHeader[ipdg]->nchmap.nbins - 1) {
          interpolatedEff = (1.5f - fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff + (-0.5f + fraction) * mLUTEntry[ipdg][inch + 1][irad][ieta][ipt]->eff;
        } else {
          interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
        }
      }
      if (mWhatEfficiency == 2) {
        if (inch < mLUTHeader[ipdg]->nchmap.nbins - 1) {
          interpolatedEff = (1.5f - fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2 + (-0.5f + fraction) * mLUTEntry[ipdg][inch + 1][irad][ieta][ipt]->eff2;
        } else {
          interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
        }
      }
    } else {
      float comparisonValue = mLUTHeader[ipdg]->nchmap.log ? log10(nch) : nch;
      if (mWhatEfficiency == 1) {
        if (inch > 0 && comparisonValue < mLUTHeader[ipdg]->nchmap.max) {
          interpolatedEff = (0.5f + fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff + (0.5f - fraction) * mLUTEntry[ipdg][inch - 1][irad][ieta][ipt]->eff;
        } else {
          interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
        }
      }
      if (mWhatEfficiency == 2) {
        if (inch > 0 && comparisonValue < mLUTHeader[ipdg]->nchmap.max) {
          interpolatedEff = (0.5f + fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2 + (0.5f - fraction) * mLUTEntry[ipdg][inch - 1][irad][ieta][ipt]->eff2;
        } else {
          interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
        }
      }
    }
  } else {
    if (mWhatEfficiency == 1)
      interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
    if (mWhatEfficiency == 2)
      interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
  }
  return mLUTEntry[ipdg][inch][irad][ieta][ipt];
} //;

/*****************************************************************/

bool TrackSmearer::smearTrack(O2Track& o2track, lutEntry_t* lutEntry, float interpolatedEff)
{
  bool isReconstructed = true;
  // generate efficiency
  if (mUseEfficiency) {
    auto eff = 0.;
    if (mWhatEfficiency == 1)
      eff = lutEntry->eff;
    if (mWhatEfficiency == 2)
      eff = lutEntry->eff2;
    if (mInterpolateEfficiency)
      eff = interpolatedEff;
    if (gRandom->Uniform() > eff)
      isReconstructed = false;
  }

  // return false already now in case not reco'ed
  if (!isReconstructed && mSkipUnreconstructed)
    return false;

  // transform params vector and smear
  double params_[5];
  for (int i = 0; i < 5; ++i) {
    double val = 0.;
    for (int j = 0; j < 5; ++j)
      val += lutEntry->eigvec[j][i] * o2track.getParam(j);
    params_[i] = gRandom->Gaus(val, sqrt(lutEntry->eigval[i]));
  }
  // transform back params vector
  for (int i = 0; i < 5; ++i) {
    double val = 0.;
    for (int j = 0; j < 5; ++j)
      val += lutEntry->eiginv[j][i] * params_[j];
    o2track.setParam(val, i);
  }
  // should make a sanity check that par[2] sin(phi) is in [-1, 1]
  if (fabs(o2track.getParam(2)) > 1.) {
    std::cout << " --- smearTrack failed sin(phi) sanity check: " << o2track.getParam(2) << std::endl;
  }
  // set covariance matrix
  for (int i = 0; i < 15; ++i)
    o2track.setCov(lutEntry->covm[i], i);
  return isReconstructed;
}

/*****************************************************************/

bool TrackSmearer::smearTrack(O2Track& o2track, int pdg, float nch)
{

  auto pt = o2track.getPt();
  if (abs(pdg) == 1000020030) {
    pt *= 2.f;
  }
  auto eta = o2track.getEta();
  float interpolatedEff = 0.0f;
  auto lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, interpolatedEff);
  if (!lutEntry || !lutEntry->valid)
    return false;
  return smearTrack(o2track, lutEntry, interpolatedEff);
}

/*****************************************************************/
// relative uncertainty on pt
double TrackSmearer::getPtRes(int pdg, float nch, float eta, float pt)
{
  float dummy = 0.0f;
  auto lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto val = sqrt(lutEntry->covm[14]) * lutEntry->pt;
  return val;
}

/*****************************************************************/
// relative uncertainty on eta
double TrackSmearer::getEtaRes(int pdg, float nch, float eta, float pt)
{
  float dummy = 0.0f;
  auto lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto sigmatgl = sqrt(lutEntry->covm[9]);                   // sigmatgl2
  auto etaRes = fabs(sin(2.0 * atan(exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  etaRes /= lutEntry->eta;                                   // relative uncertainty
  return etaRes;
}
/*****************************************************************/
// absolute uncertainty on pt
double TrackSmearer::getAbsPtRes(int pdg, float nch, float eta, float pt)
{
  float dummy = 0.0f;
  auto lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto val = sqrt(lutEntry->covm[14]) * pow(lutEntry->pt, 2);
  return val;
}

/*****************************************************************/
// absolute uncertainty on eta
double TrackSmearer::getAbsEtaRes(int pdg, float nch, float eta, float pt)
{
  float dummy = 0.0f;
  auto lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto sigmatgl = sqrt(lutEntry->covm[9]);                   // sigmatgl2
  auto etaRes = fabs(sin(2.0 * atan(exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  return etaRes;
}
/*****************************************************************/
// efficiency
double TrackSmearer::getEfficiency(int pdg, float nch, float eta, float pt)
{
  float efficiency = 0.0f;
  getLUTEntry(pdg, nch, 0., eta, pt, efficiency);
  return efficiency;
}
/*****************************************************************/
// Only in DelphesO2
// bool TrackSmearer::smearTrack(Track& track, bool atDCA)
// {

//   O2Track o2track;
//   TrackUtils::convertTrackToO2Track(track, o2track, atDCA);
//   int pdg = track.PID;
//   float nch = mdNdEta; // use locally stored dNch/deta for the time being
//   if (!smearTrack(o2track, pdg, nch))
//     return false;
//   TrackUtils::convertO2TrackToTrack(o2track, track, atDCA);
//   return true;

// #if 0
//   auto lutEntry = getLUTEntry(track.PID, 0., 0., track.Eta, track.PT);
//   if (!lutEntry)
//     return;

//   O2Track o2track;
//   TrackUtils::convertTrackToO2Track(track, o2track, atDCA);
//   smearTrack(o2track, lutEntry);
//   TrackUtils::convertO2TrackToTrack(o2track, track, atDCA);
// #endif
// }

/*****************************************************************/

} // namespace delphes
} // namespace o2
