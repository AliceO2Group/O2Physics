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
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it
///

#include "TRandom.h"
#include <map>

using namespace o2;
using namespace o2::framework;

// Imported from DelphesO2
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

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
  };
  int find(float val)
  {
    float width = (max - min) / nbins;
    int bin;
    if (log)
      bin = (int)((log10(val) - min) / width);
    else
      bin = (int)((val - min) / width);
    if (bin < 0)
      return 0;
    if (bin > nbins - 1)
      return nbins - 1;
    return bin;
  };
  void print() { printf("nbins = %d, min = %f, max = %f, log = %s \n", nbins, min, max, log ? "on" : "off"); };
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
  };
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
  };
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
  float eigvec[5][5] = {0.};
  float eiginv[5][5] = {0.};
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

class TrackSmearer
{

 public:
  TrackSmearer() = default;
  ~TrackSmearer() = default;

  /** LUT methods **/
  bool loadTable(int pdg, const char* filename, bool forceReload = false)
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
  void useEfficiency(bool val) { mUseEfficiency = val; };
  void setWhatEfficiency(int val) { mWhatEfficiency = val; };
  lutHeader_t* getLUTHeader(int pdg) { return mLUTHeader[getIndexPDG(pdg)]; };
  lutEntry_t* getLUTEntry(int pdg, float nch, float radius, float eta, float pt)
  {
    auto ipdg = getIndexPDG(pdg);
    if (!mLUTHeader[ipdg])
      return nullptr;
    auto inch = mLUTHeader[ipdg]->nchmap.find(nch);
    auto irad = mLUTHeader[ipdg]->radmap.find(radius);
    auto ieta = mLUTHeader[ipdg]->etamap.find(eta);
    auto ipt = mLUTHeader[ipdg]->ptmap.find(pt);
    return mLUTEntry[ipdg][inch][irad][ieta][ipt];
  }

  using O2Track = o2::track::TrackParCov;
  bool smearTrack(O2Track& o2track, lutEntry_t* lutEntry)
  {
    // generate efficiency
    if (mUseEfficiency) {
      auto eff = 0.;
      if (mWhatEfficiency == 1)
        eff = lutEntry->eff;
      if (mWhatEfficiency == 2)
        eff = lutEntry->eff2;
      if (gRandom->Uniform() > eff)
        return false;
    }
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
    return true;
  }
  bool smearTrack(O2Track& o2track, int pid, float nch)
  {

    auto pt = o2track.getPt();
    if (abs(pid) == 1000020030) {
      pt *= 2.f;
    }
    auto eta = o2track.getEta();
    auto lutEntry = getLUTEntry(pid, nch, 0., eta, pt);
    if (!lutEntry || !lutEntry->valid)
      return false;
    return smearTrack(o2track, lutEntry);
  }

  int getIndexPDG(int pdg)
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
      default:
        return 2; // Default: pion
    };
  };

  void setdNdEta(float val) { mdNdEta = val; };

 protected:
  static constexpr unsigned int nLUTs = 8; // Number of LUT available
  lutHeader_t* mLUTHeader[nLUTs] = {nullptr};
  lutEntry_t***** mLUTEntry[nLUTs] = {nullptr};
  bool mUseEfficiency = true;
  int mWhatEfficiency = 1;
  float mdNdEta = 1600.;
};
