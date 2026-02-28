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
/// \file DelphesO2TrackSmearer.cxx
/// \author Roberto Preghenella
/// \brief Porting to O2Physics of DelphesO2 code.
///        Relevant sources:
///                 DelphesO2/src/lutCovm.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutCovm.hh
///                 DelphesO2/src/TrackSmearer.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.cc
///                 DelphesO2/src/TrackSmearer.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.hh
/// @email: preghenella@bo.infn.it
///

#include "ALICE3/Core/DelphesO2TrackSmearer.h"

#include "ALICE3/Core/GeometryContainer.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <map>
#include <string>

namespace o2
{
namespace delphes
{

/*****************************************************************/

int DelphesO2TrackSmearer::getIndexPDG(const int pdg)
{
  switch (abs(pdg)) {
    case 11:
      return 0; // Electron
    case 13:
      return 1; // Muon
    case 211:
      return 2; // Pion
    case 321:
      return 3; // Kaon
    case 2212:
      return 4; // Proton
    case 1000010020:
      return 5; // Deuteron
    case 1000010030:
      return 6; // Triton
    case 1000020030:
      return 7; // Helium3
    case 1000020040:
      return 8; // Alphas
    default:
      return 2; // Default: pion
  }
}

const char* DelphesO2TrackSmearer::getParticleName(int pdg)
{
  switch (abs(pdg)) {
    case 11:
      return "electron";
    case 13:
      return "muon";
    case 211:
      return "pion";
    case 321:
      return "kaon";
    case 2212:
      return "proton";
    case 1000010020:
      return "deuteron";
    case 1000010030:
      return "triton";
    case 1000020030:
      return "helium3";
    case 1000020040:
      return "alpha";
    default:
      return "pion"; // Default: pion
  }
}

bool DelphesO2TrackSmearer::loadTable(int pdg, const char* filename, bool forceReload)
{
  if (!filename || filename[0] == '\0') {
    LOG(info) << " --- No LUT file provided for PDG " << pdg << ". Skipping load.";
    return false;
  }
  const auto ipdg = getIndexPDG(pdg);
  LOGF(info, "Will load %s lut file ..: '%s'", getParticleName(pdg), filename);
  if (mLUTHeader[ipdg] && !forceReload) {
    LOG(info) << " --- LUT table for PDG " << pdg << " has been already loaded with index " << ipdg << std::endl;
    return false;
  }

  const std::string localFilename = o2::fastsim::GeometryEntry::accessFile(filename, "./.ALICE3/LUTs/", mCcdbManager, 10);
  mLUTHeader[ipdg] = new lutHeader_t;

  std::ifstream lutFile(localFilename, std::ifstream::binary);
  if (!lutFile.is_open()) {
    LOG(info) << " --- cannot open covariance matrix file for PDG " << pdg << ": " << localFilename << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  lutFile.read(reinterpret_cast<char*>(mLUTHeader[ipdg]), sizeof(lutHeader_t));
  if (lutFile.gcount() != sizeof(lutHeader_t)) {
    LOG(info) << " --- troubles reading covariance matrix header for PDG " << pdg << ": " << filename << std::endl;
    LOG(info) << " --- expected/detected " << sizeof(lutHeader_t) << "/" << lutFile.gcount() << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  if (mLUTHeader[ipdg]->version != LUTCOVM_VERSION) {
    LOG(info) << " --- LUT header version mismatch: expected/detected = " << LUTCOVM_VERSION << "/" << mLUTHeader[ipdg]->version << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  bool specialPdgCase = false;
  switch (pdg) {                         // Handle special cases
    case o2::constants::physics::kAlpha: // Special case: Allow Alpha particles to use He3 LUT
      specialPdgCase = (mLUTHeader[ipdg]->pdg == o2::constants::physics::kHelium3);
      if (specialPdgCase)
        LOG(info)
          << " --- Alpha particles (PDG " << pdg << ") will use He3 LUT data (PDG " << mLUTHeader[ipdg]->pdg << ")" << std::endl;
      break;
    default:
      break;
  }
  if (mLUTHeader[ipdg]->pdg != pdg && !specialPdgCase) {
    LOG(info) << " --- LUT header PDG mismatch: expected/detected = " << pdg << "/" << mLUTHeader[ipdg]->pdg << std::endl;
    delete mLUTHeader[ipdg];
    mLUTHeader[ipdg] = nullptr;
    return false;
  }
  const int nnch = mLUTHeader[ipdg]->nchmap.nbins;
  const int nrad = mLUTHeader[ipdg]->radmap.nbins;
  const int neta = mLUTHeader[ipdg]->etamap.nbins;
  const int npt = mLUTHeader[ipdg]->ptmap.nbins;
  mLUTEntry[ipdg] = new DelphesO2TrackSmearer::lutEntry_t****[nnch];
  for (int inch = 0; inch < nnch; ++inch) {
    mLUTEntry[ipdg][inch] = new DelphesO2TrackSmearer::lutEntry_t***[nrad];
    for (int irad = 0; irad < nrad; ++irad) {
      mLUTEntry[ipdg][inch][irad] = new DelphesO2TrackSmearer::lutEntry_t**[neta];
      for (int ieta = 0; ieta < neta; ++ieta) {
        mLUTEntry[ipdg][inch][irad][ieta] = new DelphesO2TrackSmearer::lutEntry_t*[npt];
        for (int ipt = 0; ipt < npt; ++ipt) {
          mLUTEntry[ipdg][inch][irad][ieta][ipt] = new DelphesO2TrackSmearer::lutEntry_t;
          lutFile.read(reinterpret_cast<char*>(mLUTEntry[ipdg][inch][irad][ieta][ipt]), sizeof(DelphesO2TrackSmearer::lutEntry_t));
          if (lutFile.gcount() != sizeof(DelphesO2TrackSmearer::lutEntry_t)) {
            LOG(info) << " --- troubles reading covariance matrix entry for PDG " << pdg << ": " << localFilename << std::endl;
            LOG(info) << " --- expected/detected " << sizeof(lutHeader_t) << "/" << lutFile.gcount() << std::endl;
            return false;
          }
        }
      }
    }
  }
  LOG(info) << " --- read covariance matrix table for PDG " << pdg << ": " << filename << std::endl;
  mLUTHeader[ipdg]->print();

  lutFile.close();
  return true;
}

/*****************************************************************/

DelphesO2TrackSmearer::lutEntry_t* DelphesO2TrackSmearer::getLUTEntry(const int pdg, const float nch, const float radius, const float eta, const float pt, float& interpolatedEff)
{
  const int ipdg = getIndexPDG(pdg);
  if (!mLUTHeader[ipdg]) {
    return nullptr;
  }

  auto inch = mLUTHeader[ipdg]->nchmap.find(nch);
  auto irad = mLUTHeader[ipdg]->radmap.find(radius);
  auto ieta = mLUTHeader[ipdg]->etamap.find(eta);
  auto ipt = mLUTHeader[ipdg]->ptmap.find(pt);

  // Interpolate if requested
  auto fraction = mLUTHeader[ipdg]->nchmap.fracPositionWithinBin(nch);
  if (mInterpolateEfficiency) {
    static constexpr float kFractionThreshold = 0.5f;
    if (fraction > kFractionThreshold) {
      switch (mWhatEfficiency) {
        case 1:
          if (inch < mLUTHeader[ipdg]->nchmap.nbins - 1) {
            interpolatedEff = (1.5f - fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff + (-0.5f + fraction) * mLUTEntry[ipdg][inch + 1][irad][ieta][ipt]->eff;
          } else {
            interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
          }
          break;
        case 2:
          if (inch < mLUTHeader[ipdg]->nchmap.nbins - 1) {
            interpolatedEff = (1.5f - fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2 + (-0.5f + fraction) * mLUTEntry[ipdg][inch + 1][irad][ieta][ipt]->eff2;
          } else {
            interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
          }
          break;
        default:
          LOG(fatal) << " --- getLUTEntry: unknown efficiency type " << mWhatEfficiency;
      }
    } else {
      float comparisonValue = mLUTHeader[ipdg]->nchmap.log ? std::log10(nch) : nch;
      switch (mWhatEfficiency) {
        case 1:
          if (inch > 0 && comparisonValue < mLUTHeader[ipdg]->nchmap.max) {
            interpolatedEff = (0.5f + fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff + (0.5f - fraction) * mLUTEntry[ipdg][inch - 1][irad][ieta][ipt]->eff;
          } else {
            interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
          }
          break;
        case 2:
          if (inch > 0 && comparisonValue < mLUTHeader[ipdg]->nchmap.max) {
            interpolatedEff = (0.5f + fraction) * mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2 + (0.5f - fraction) * mLUTEntry[ipdg][inch - 1][irad][ieta][ipt]->eff2;
          } else {
            interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
          }
          break;
        default:
          LOG(fatal) << " --- getLUTEntry: unknown efficiency type " << mWhatEfficiency;
      }
    }
  } else {
    switch (mWhatEfficiency) {
      case 1:
        interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff;
        break;
      case 2:
        interpolatedEff = mLUTEntry[ipdg][inch][irad][ieta][ipt]->eff2;
        break;
      default:
        LOG(fatal) << " --- getLUTEntry: unknown efficiency type " << mWhatEfficiency;
    }
  }
  return mLUTEntry[ipdg][inch][irad][ieta][ipt];
} //;

/*****************************************************************/

bool DelphesO2TrackSmearer::smearTrack(o2::track::TrackParCov& o2track, DelphesO2TrackSmearer::lutEntry_t* lutEntry, float interpolatedEff)
{
  bool isReconstructed = true;
  // generate efficiency
  if (mUseEfficiency) {
    auto eff = 0.;
    switch (mWhatEfficiency) {
      case 1:
        eff = lutEntry->eff;
        break;
      case 2:
        eff = lutEntry->eff2;
        break;
    }
    if (mInterpolateEfficiency)
      eff = interpolatedEff;
    if (gRandom->Uniform() > eff)
      isReconstructed = false;
  }

  // return false already now in case not reco'ed
  if (!isReconstructed && mSkipUnreconstructed)
    return false;

  // transform params vector and smear
  static constexpr int kParSize = 5;
  double params[kParSize];
  for (int i = 0; i < kParSize; ++i) {
    double val = 0.;
    for (int j = 0; j < kParSize; ++j)
      val += lutEntry->eigvec[j][i] * o2track.getParam(j);
    params[i] = gRandom->Gaus(val, std::sqrt(lutEntry->eigval[i]));
  }
  // transform back params vector
  for (int i = 0; i < kParSize; ++i) {
    double val = 0.;
    for (int j = 0; j < kParSize; ++j)
      val += lutEntry->eiginv[j][i] * params[j];
    o2track.setParam(val, i);
  }
  // should make a sanity check that par[2] sin(phi) is in [-1, 1]
  if (std::fabs(o2track.getParam(2)) > 1.) {
    LOG(info) << " --- smearTrack failed sin(phi) sanity check: " << o2track.getParam(2) << std::endl;
  }
  // set covariance matrix
  static constexpr int kCovMatSize = 15;
  for (int i = 0; i < kCovMatSize; ++i)
    o2track.setCov(lutEntry->covm[i], i);
  return isReconstructed;
}

/*****************************************************************/

bool DelphesO2TrackSmearer::smearTrack(o2::track::TrackParCov& o2track, int pdg, float nch)
{
  auto pt = o2track.getPt();
  switch (pdg) {
    case o2::constants::physics::kHelium3:
    case -o2::constants::physics::kHelium3:
      pt *= 2.f;
      break;
  }
  auto eta = o2track.getEta();
  float interpolatedEff = 0.0f;
  DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, interpolatedEff);
  if (!lutEntry || !lutEntry->valid)
    return false;
  return smearTrack(o2track, lutEntry, interpolatedEff);
}

/*****************************************************************/
// relative uncertainty on pt
double DelphesO2TrackSmearer::getPtRes(const int pdg, const float nch, const float eta, const float pt)
{
  float dummy = 0.0f;
  DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt;
  return val;
}

/*****************************************************************/
// relative uncertainty on eta
double DelphesO2TrackSmearer::getEtaRes(const int pdg, const float nch, const float eta, const float pt)
{
  float dummy = 0.0f;
  DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto sigmatgl = std::sqrt(lutEntry->covm[9]);                                  // sigmatgl2
  auto etaRes = std::fabs(std::sin(2.0 * std::atan(std::exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  etaRes /= lutEntry->eta;                                                       // relative uncertainty
  return etaRes;
}
/*****************************************************************/
// absolute uncertainty on pt
double DelphesO2TrackSmearer::getAbsPtRes(const int pdg, const float nch, const float eta, const float pt)
{
  float dummy = 0.0f;
  DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto val = std::sqrt(lutEntry->covm[14]) * lutEntry->pt * lutEntry->pt;
  return val;
}

/*****************************************************************/
// absolute uncertainty on eta
double DelphesO2TrackSmearer::getAbsEtaRes(const int pdg, const float nch, const float eta, const float pt)
{
  float dummy = 0.0f;
  DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(pdg, nch, 0., eta, pt, dummy);
  auto sigmatgl = std::sqrt(lutEntry->covm[9]);                                  // sigmatgl2
  auto etaRes = std::fabs(std::sin(2.0 * std::atan(std::exp(-eta)))) * sigmatgl; // propagate tgl to eta uncertainty
  return etaRes;
}
/*****************************************************************/
// efficiency
double DelphesO2TrackSmearer::getEfficiency(const int pdg, const float nch, const float eta, const float pt)
{
  float efficiency = 0.0f;
  getLUTEntry(pdg, nch, 0., eta, pt, efficiency);
  return efficiency;
}
/*****************************************************************/
// Only in DelphesO2
// bool DelphesO2TrackSmearer::smearTrack(Track& track, bool atDCA)
// {

//   o2::track::TrackParCov o2track;
//   TrackUtils::convertTrackToO2Track(track, o2track, atDCA);
//   int pdg = track.PID;
//   float nch = mdNdEta; // use locally stored dNch/deta for the time being
//   if (!smearTrack(o2track, pdg, nch))
//     return false;
//   TrackUtils::convertO2TrackToTrack(o2track, track, atDCA);
//   return true;

// #if 0
//   DelphesO2TrackSmearer::lutEntry_t* lutEntry = getLUTEntry(track.PID, 0., 0., track.Eta, track.PT);
//   if (!lutEntry)
//     return;

//   o2::track::TrackParCov o2track;
//   TrackUtils::convertTrackToO2Track(track, o2track, atDCA);
//   smearTrack(o2track, lutEntry);
//   TrackUtils::convertO2TrackToTrack(o2track, track, atDCA);
// #endif
// }

/*****************************************************************/

} // namespace delphes
} // namespace o2
