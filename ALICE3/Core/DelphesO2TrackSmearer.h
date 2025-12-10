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
/// @file DelphesO2TrackSmearer.h
/// @brief Porting to O2Physics of DelphesO2 code.
///        Minimal changes have been made to the original code for adaptation purposes, formatting and commented parts have been considered.
///        Relevant sources:
///                 DelphesO2/src/lutCovm.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/lutCovm.hh
///                 DelphesO2/src/TrackSmearer.cc https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.cc
///                 DelphesO2/src/TrackSmearer.hh https://github.com/AliceO2Group/DelphesO2/blob/master/src/TrackSmearer.hh
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it
///

#ifndef ALICE3_CORE_DELPHESO2TRACKSMEARER_H_
#define ALICE3_CORE_DELPHESO2TRACKSMEARER_H_

#include <CCDB/BasicCCDBManager.h>
#include <ReconstructionDataFormats/Track.h>

#include <TRandom.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>

///////////////////////////////
/// DelphesO2/src/lutCovm.hh //
///////////////////////////////

/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

// #pragma // once
#define LUTCOVM_VERSION 20210801

struct map_t {
  int nbins = 1;
  float min = 0.;
  float max = 1.e6;
  bool log = false;
  float eval(int bin)
  {
    float width = (max - min) / nbins;
    float val = min + (bin + 0.5) * width;
    if (log)
      return pow(10., val);
    return val;
  }
  // function needed to interpolate some dimensions
  float fracPositionWithinBin(float val)
  {
    float width = (max - min) / nbins;
    int bin;
    float returnVal = 0.5f;
    if (log) {
      bin = static_cast<int>((log10(val) - min) / width);
      returnVal = ((log10(val) - min) / width) - bin;
    } else {
      bin = static_cast<int>((val - min) / width);
      returnVal = val / width - bin;
    }
    return returnVal;
  }

  int find(float val)
  {
    float width = (max - min) / nbins;
    int bin;
    if (log)
      bin = static_cast<int>((log10(val) - min) / width); // Changed due to MegaLinter error.
    // bin = (int)((log10(val) - min) / width); // Original line.
    else
      bin = static_cast<int>((val - min) / width); // Changed due to MegaLinter error.
    // bin = (int)((val - min) / width); // Original line.
    if (bin < 0)
      return 0;
    if (bin > nbins - 1)
      return nbins - 1;
    return bin;
  } //;
  void print() { printf("nbins = %d, min = %f, max = %f, log = %s \n", nbins, min, max, log ? "on" : "off"); } //;
};

struct lutHeader_t {
  int version = LUTCOVM_VERSION;
  int pdg = 0;
  float mass = 0.;
  float field = 0.;
  map_t nchmap;
  map_t radmap;
  map_t etamap;
  map_t ptmap;
  bool check_version()
  {
    return (version == LUTCOVM_VERSION);
  } //;
  void print()
  {
    printf(" version: %d \n", version);
    printf("     pdg: %d \n", pdg);
    printf("   field: %f \n", field);
    printf("  nchmap: ");
    nchmap.print();
    printf("  radmap: ");
    radmap.print();
    printf("  etamap: ");
    etamap.print();
    printf("   ptmap: ");
    ptmap.print();
  } //;
};

struct lutEntry_t {
  float nch = 0.;
  float eta = 0.;
  float pt = 0.;
  bool valid = false;
  float eff = 0.;
  float eff2 = 0.;
  float itof = 0.;
  float otof = 0.;
  float covm[15] = {0.};
  float eigval[5] = {0.};
  float eigvec[5][5] = {{0.}};
  float eiginv[5][5] = {{0.}};
  void print()
  {
    printf(" --- lutEntry: pt = %f, eta = %f (%s)\n", pt, eta, valid ? "valid" : "not valid");
    printf("     efficiency: %f\n", eff);
    printf("     covMatix: ");
    int k = 0;
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < i + 1; ++j)
        printf("% e ", covm[k++]);
      printf("\n               ");
    }
    printf("\n");
  }
};

////////////////////////////////////
/// DelphesO2/src/TrackSmearer.hh //
////////////////////////////////////

/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

// #ifndef _DelphesO2_TrackSmearer_h_
// #define _DelphesO2_TrackSmearer_h_

// #include "ReconstructionDataFormats/Track.h"
// #include "classes/DelphesClasses.h"
// #include "lutCovm.hh"
// #include <map>

using O2Track = o2::track::TrackParCov;

namespace o2
{
namespace delphes
{

class TrackSmearer
{

 public:
  TrackSmearer() = default;
  ~TrackSmearer() = default;

  /** LUT methods **/
  bool loadTable(int pdg, const char* filename, bool forceReload = false);
  bool hasTable(int pdg) { return (mLUTHeader[getIndexPDG(pdg)] != nullptr); } //;
  void useEfficiency(bool val) { mUseEfficiency = val; }                       //;
  void interpolateEfficiency(bool val) { mInterpolateEfficiency = val; }       //;
  void skipUnreconstructed(bool val) { mSkipUnreconstructed = val; }           //;
  void setWhatEfficiency(int val) { mWhatEfficiency = val; }                   //;
  lutHeader_t* getLUTHeader(int pdg) { return mLUTHeader[getIndexPDG(pdg)]; }  //;
  lutEntry_t* getLUTEntry(const int pdg, const float nch, const float radius, const float eta, const float pt, float& interpolatedEff);

  bool smearTrack(O2Track& o2track, lutEntry_t* lutEntry, float interpolatedEff);
  bool smearTrack(O2Track& o2track, int pdg, float nch);
  // bool smearTrack(Track& track, bool atDCA = true); // Only in DelphesO2
  double getPtRes(const int pdg, const float nch, const float eta, const float pt);
  double getEtaRes(const int pdg, const float nch, const float eta, const float pt);
  double getAbsPtRes(const int pdg, const float nch, const float eta, const float pt);
  double getAbsEtaRes(const int pdg, const float nch, const float eta, const float pt);
  double getEfficiency(const int pdg, const float nch, const float eta, const float pt);

  int getIndexPDG(const int pdg)
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

  const char* getParticleName(int pdg)
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
  void setdNdEta(float val) { mdNdEta = val; }                                 //;
  void setCcdbManager(o2::ccdb::BasicCCDBManager* mgr) { mCcdbManager = mgr; } //;
  void setCleanupDownloadedFile(bool val) { mCleanupDownloadedFile = val; }    //;

 protected:
  static constexpr unsigned int nLUTs = 9; // Number of LUT available
  lutHeader_t* mLUTHeader[nLUTs] = {nullptr};
  lutEntry_t***** mLUTEntry[nLUTs] = {nullptr};
  bool mUseEfficiency = true;
  bool mInterpolateEfficiency = false;
  bool mSkipUnreconstructed = true; // don't smear tracks that are not reco'ed
  int mWhatEfficiency = 1;
  float mdNdEta = 1600.;

 private:
  o2::ccdb::BasicCCDBManager* mCcdbManager = nullptr;
  bool mCleanupDownloadedFile = true;
};

} // namespace delphes
} // namespace o2

// #endif /** _DelphesO2_TrackSmearer_h_ **/

namespace o2::delphes
{
using DelphesO2TrackSmearer = TrackSmearer;
}
#endif // ALICE3_CORE_DELPHESO2TRACKSMEARER_H_
