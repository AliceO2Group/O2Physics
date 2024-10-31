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

#ifndef PWGLF_UTILS_V0SELECTIONGROUP_H_
#define PWGLF_UTILS_V0SELECTIONGROUP_H_

#include <iosfwd>
#include <Rtypes.h>
#include <TMath.h>
#include "v0SelectionBits.h"
#include "Framework/Logger.h"
#include "CommonConstants/PhysicsConstants.h"

class v0SelectionGroup
{
 public:
  v0SelectionGroup()
    : rapidityCut{0.5f}, daughterEtaCut{0.8f}, v0cospa{0.97f}, dcav0dau{1.0f}, dcanegtopv{0.05f}, dcapostopv{0.05f}, v0radius{1.2f}, v0radiusMax{1e+5}, minTPCrows{70}, minITSclusters{-1}, skipTPConly{false}, requirePosITSonly{false}, requireNegITSonly{false}, TpcPidNsigmaCut{5.0f}, TofPidNsigmaCutLaPr{1e+6}, TofPidNsigmaCutLaPi{1e+6}, TofPidNsigmaCutK0Pi{1e+6}, maxDeltaTimeProton{1e+9}, maxDeltaTimePion{1e+9}, lifetimeCutK0Short{20.0f}, lifetimeCutLambda{20.0f}, armPodCut{5.0f}
  {
    // constructor
  }

  void provideMasks(uint64_t& maskTopological, uint64_t& maskTrackProperties, uint64_t& maskK0ShortSpecific, uint64_t& maskLambdaSpecific, uint64_t& maskAntiLambdaSpecific) const;
  bool verifyMask(uint64_t bitmap, uint64_t mask) const;

  float getRapidityCut() const { return rapidityCut; }
  float getDaughterEtaCut() const { return daughterEtaCut; }

  float getv0cospa() const { return v0cospa; }
  float getdcav0dau() const { return dcav0dau; }
  float getdcanegtopv() const { return dcanegtopv; }
  float getdcapostopv() const { return dcapostopv; }
  float getv0radius() const { return v0radius; }
  float getv0radiusMax() const { return v0radiusMax; }

  int getminTPCrows() const { return minTPCrows; }
  int getminITSclusters() const { return minITSclusters; }
  bool getskipTPConly() const { return skipTPConly; }
  bool getrequirePosITSonly() const { return requirePosITSonly; }
  bool getrequireNegITSonly() const { return requireNegITSonly; }

  float getTpcPidNsigmaCut() const { return TpcPidNsigmaCut; }
  float getTofPidNsigmaCutLaPr() const { return TofPidNsigmaCutLaPr; }
  float getTofPidNsigmaCutLaPi() const { return TofPidNsigmaCutLaPi; }
  float getTofPidNsigmaCutK0Pi() const { return TofPidNsigmaCutK0Pi; }

  float getmaxDeltaTimeProton() const { return maxDeltaTimeProton; }
  float getmaxDeltaTimePion() const { return maxDeltaTimePion; }

  float getlifetimeCutK0Short() const { return lifetimeCutK0Short; }
  float getlifetimeCutLambda() const { return lifetimeCutLambda; }

  float getarmPodCut() const { return armPodCut; }

  // Helper to print out selections
  void PrintSelections() const;

 private:
  // Phase space
  float rapidityCut;
  float daughterEtaCut;
  // Topology
  float v0cospa;
  float dcav0dau;
  float dcanegtopv;
  float dcapostopv;
  float v0radius;
  float v0radiusMax;

  // Track quality
  int minTPCrows;
  int minITSclusters;
  bool skipTPConly;
  bool requirePosITSonly;
  bool requireNegITSonly;

  // PID
  float TpcPidNsigmaCut;
  float TofPidNsigmaCutLaPr;
  float TofPidNsigmaCutLaPi;
  float TofPidNsigmaCutK0Pi;

  // PID precursor (compatibility only);
  float maxDeltaTimeProton;
  float maxDeltaTimePion;

  // Lifetime
  float lifetimeCutK0Short;
  float lifetimeCutLambda;

  // Armenteros
  float armPodCut;

  ClassDefNV(v0SelectionGroup, 1);
};

#endif // PWGLF_UTILS_V0SELECTIONGROUP_H_
