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

#include "EMCPhotonCut.h"

#include "PWGJE/DataModel/EMCALClusters.h"

#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <Rtypes.h>

#include <string>

ClassImp(EMCPhotonCut);

const char* EMCPhotonCut::mCutNames[static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts)] = {"Definition", "Energy", "NCell", "M02", "Timing", "TrackMatching", "SecTrackMatching", "Exotic"};

void EMCPhotonCut::addQAHistograms(o2::framework::HistogramRegistry* fRegistry) const
{
  if (mDoQA && fRegistry != nullptr) {
    const o2::framework::AxisSpec thAxisClusterEnergy{500, 0, 50, "#it{E}_{cls} (GeV)"};
    const o2::framework::AxisSpec thAxisMomentum{250, 0., 25., "#it{p}_{T} (GeV/#it{c})"};
    const o2::framework::AxisSpec thAxisDEta{200, -0.1, 0.1, "#Delta#eta"};
    const o2::framework::AxisSpec thAxisDPhi{200, -0.1, 0.1, "#Delta#varphi (rad)"};
    const o2::framework::AxisSpec thAxisEnergy{500, 0., 50., "#it{E} (GeV)"};
    const o2::framework::AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
    const o2::framework::AxisSpec thAxisPhi{500, 0, o2::constants::math::TwoPI, "#varphi (rad)"};
    const o2::framework::AxisSpec thAxisNCell{51, -0.5, 50.5, "#it{N}_{cell}"};
    const o2::framework::AxisSpec thAxisM02{200, 0, 2.0, "#it{M}_{02}"};
    const o2::framework::AxisSpec thAxisTime{300, -150, +150, "#it{t}_{cls} (ns)"};
    const o2::framework::AxisSpec thAxisEoverP{400, 0, 10., "#it{E}_{cls}/#it{p}_{track} (#it{c})"};

    fRegistry->add("QA/Cluster/before/hE", "E_{cluster};#it{E}_{cluster} (GeV);#it{N}_{cluster}", o2::framework::HistType::kTH1D, {thAxisClusterEnergy}, true);
    fRegistry->add("QA/Cluster/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{cluster}", o2::framework::HistType::kTH1D, {thAxisClusterEnergy}, true);
    fRegistry->add("QA/Cluster/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", o2::framework::HistType::kTH1D, {{1001, -0.5f, 1000.5f}}, true);
    fRegistry->add("QA/Cluster/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisEta, thAxisPhi}, true);
    fRegistry->add("QA/Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{E}_{cluster} (GeV)", o2::framework::HistType::kTH2F, {thAxisNCell, thAxisClusterEnergy}, true);
    fRegistry->add("QA/Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{E}_{cluster} (GeV)", o2::framework::HistType::kTH2F, {thAxisM02, thAxisClusterEnergy}, true);
    fRegistry->add("QA/Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{E}_{cluster} (GeV)", o2::framework::HistType::kTH2F, {thAxisTime, thAxisClusterEnergy}, true);

    fRegistry->addClone("QA/Cluster/before/", "QA/Cluster/after/");

    auto hClusterQualityCuts = fRegistry->add<TH2>("QA/Cluster/hClusterQualityCuts", "Energy at which clusters are removed by a given cut;;#it{E} (GeV)", o2::framework::HistType::kTH2F, {{static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 2, -0.5, static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 1.5}, thAxisClusterEnergy}, true);
    hClusterQualityCuts->GetXaxis()->SetBinLabel(1, "In");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(2, "Definition");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(3, "Energy");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(4, "NCell");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(5, "M02");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(6, "Timing");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(7, "TM");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(8, "Sec. TM");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(9, "Exotic");
    hClusterQualityCuts->GetXaxis()->SetBinLabel(10, "Out");

    fRegistry->add("QA/Cluster/hTrackdEtadPhi", "d#eta vs. d#varphi of matched tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDEta, thAxisDPhi}, true);
    fRegistry->add("QA/Cluster/hTrackdEtaPt", "d#eta vs. track pT of matched tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDEta, thAxisMomentum}, true);
    fRegistry->add("QA/Cluster/hTrackdPhiPt", "d#varphi vs. track pT of matched tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDPhi, thAxisMomentum}, true);
    fRegistry->add("QA/Cluster/hSecTrackdEtadPhi", "d#eta vs. d#varphi of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDEta, thAxisDPhi}, true);
    fRegistry->add("QA/Cluster/hSecTrackdEtaPt", "d#eta vs. track pT of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDEta, thAxisMomentum}, true);
    fRegistry->add("QA/Cluster/hSecTrackdPhiPt", "d#varphi vs. track pT of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::HistType::kTH2F, {thAxisDPhi, thAxisMomentum}, true);
  }
}

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
