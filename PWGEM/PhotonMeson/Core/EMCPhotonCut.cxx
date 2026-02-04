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

/// \file EMCPhotonCut.cxx
/// \brief source of class for emcal photon selection.
/// \author M. Hemmer, marvin.hemmer@cern.ch; N. Strangmann, nicolas.strangmann@cern.ch

#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
//
#include "PWGJE/DataModel/EMCALClusters.h"

#include <Framework/Logger.h>

#include <Rtypes.h>

#include <string>

ClassImp(EMCPhotonCut);

const char* EMCPhotonCut::mCutNames[static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts)] = {"Definition", "Energy", "NCell", "M02", "Timing", "TrackMatching", "SecTrackMatching", "Exotic"};

void EMCPhotonCut::SetClusterizer(std::string clusterDefinitionString)
{
  mDefinition = static_cast<int>(o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionString));
  LOG(info) << "EMCal Photon Cut, set cluster definition to: " << mDefinition << " (" << clusterDefinitionString << ")";
}

void EMCPhotonCut::SetMinE(float min)
{
  mMinE = min;
  LOG(info) << "EMCal Photon Cut, set minimum cluster energy: " << mMinE;
}

void EMCPhotonCut::SetMinNCell(int min)
{
  mMinNCell = min;
  LOG(info) << "EMCal Photon Cut, set minimum number of cells per cluster: " << mMinNCell;
}

void EMCPhotonCut::SetM02Range(float min, float max)
{
  mMinM02 = min;
  mMaxM02 = max;
  LOG(info) << "EMCal Photon Cut, set minimum and maximum M02: " << mMinM02 << " <= M02 <= " << mMaxM02;
}

void EMCPhotonCut::SetTimeRange(float min, float max)
{
  mMinTime = min;
  mMaxTime = max;
  LOG(info) << "EMCal Photon Cut, set cluster time range in ns: " << mMinTime << " <= t <= " << mMaxTime;
}

void EMCPhotonCut::SetMinEoverP(float min)
{
  mMinEoverP = min;
}

void EMCPhotonCut::SetUseExoticCut(bool flag)
{
  mUseExoticCut = flag;
  LOG(info) << "EMCal Photon Cut, set usage of exotic cluster cut to: " << mUseExoticCut;
}

void EMCPhotonCut::SetUseTM(bool flag)
{
  mUseTM = flag;
  LOG(info) << "EM Photon Cluster Cut, using TM cut is set to : " << mUseTM;
}

void EMCPhotonCut::SetUseSecondaryTM(bool flag)
{
  mUseSecondaryTM = flag;
  LOG(info) << "EM Photon Cluster Cut, using secondary TM cut is set to : " << mUseTM;
}

void EMCPhotonCut::SetDoQA(bool flag)
{
  mDoQA = flag;
  LOG(info) << "EM Photon Cluster Cut, QA is set to: " << mUseTM;
}

void EMCPhotonCut::print() const
{
  LOG(info) << "EMCal Photon Cut:";
  for (int i = 0; i < static_cast<int>(EMCPhotonCuts::kNCuts); i++) {
    switch (static_cast<EMCPhotonCuts>(i)) {
      case EMCPhotonCuts::kDefinition:
        LOG(info) << mCutNames[i] << " > " << mDefinition;
        break;
      case EMCPhotonCuts::kEnergy:
        LOG(info) << mCutNames[i] << " > " << mMinE;
        break;
      case EMCPhotonCuts::kNCell:
        LOG(info) << mCutNames[i] << " > " << mMinNCell;
        break;
      case EMCPhotonCuts::kM02:
        LOG(info) << mCutNames[i] << " in [" << mMinM02 << ", " << mMaxM02 << "]";
        break;
      case EMCPhotonCuts::kTiming:
        LOG(info) << mCutNames[i] << " in [" << mMinTime << ", " << mMaxTime << "]";
        break;
      // TODO: find a nice way to print TM cuts
      // case EMCPhotonCuts::kTM:
      //   LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsOverFindableClustersTPC;
      //   break;
      case EMCPhotonCuts::kExotic:
        LOG(info) << mCutNames[i] << " set to " << mUseExoticCut;
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
