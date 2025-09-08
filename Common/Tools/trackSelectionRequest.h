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

// This class serves to aggregate track selection criteria in a standardized
// configurable that can be queried at analysis time by any service that
// could pre-apply some of these selections. It also serves as an access point
// for these cuts in a local task, where these selections are also presumed
// to be applied for safety whenever operating with multiple analyses.
//
// Following a 'request-and-reply' logic, in a certain analysis topology, all
// requests from analysis tasks should be queried by service tasks and
// the actual selection should be based on a logical || of multiple requests.
// Because of this, it is particularly important that the cuts in this object
// in an analysis!

#ifndef COMMON_TOOLS_TRACKSELECTIONREQUEST_H_
#define COMMON_TOOLS_TRACKSELECTIONREQUEST_H_

#include <Rtypes.h>

#include <cmath>
#include <iosfwd>

class trackSelectionRequest
{
 public:
  trackSelectionRequest()
    : trackPhysicsType{0}, minPt{0.0}, maxPt{1e+6}, minEta{-100}, maxEta{+100}, maxDCAz{1e+6}, maxDCAxyPtDep{1e+6}, requireTPC{false}, minTPCclusters{-1}, minTPCcrossedrows{-1}, minTPCcrossedrowsoverfindable{0.0}, maxTPCFractionSharedCls{0.0}, requireITS{false}, minITSclusters{-1}, maxITSChi2percluster{1e+6}
  {
    // constructor
  }
  void setTrackPhysicsType(int trackPhysicsType_);
  int getTrackPhysicsType() const;
  void setMinPt(float minPt_);
  int getMinPt() const;
  void setMaxPt(float maxPt_);
  int getMaxPt() const;
  void setMinEta(float minEta_);
  int getMinEta() const;
  void setMaxEta(float maxEta_);
  int getMaxEta() const;

  void setMaxDCAz(float maxDCAz_);
  int getMaxDCAz() const;
  void setMaxDCAxyPtDep(float maxDCAxyPtDep_);
  int getMaxDCAxyPtDep() const;

  void setRequireTPC(bool requireTPC_);
  bool getRequireTPC() const;
  void setMinTPCClusters(int minTPCclusters_);
  int getMinTPCClusters() const;
  void setMinTPCCrossedRows(int minTPCCrossedRows_);
  int getMinTPCCrossedRows() const;
  void setMinTPCCrossedRowsOverFindable(float minTPCCrossedRowsOverFindable_);
  int getMinTPCCrossedRowsOverFindable() const;
  void setMaxTPCFractionSharedCls(float maxTPCFractionSharedCls_);
  int getMaxTPCFractionSharedCls() const;

  void setRequireITS(bool requireITS_);
  bool getRequireITS() const;
  void setMinITSClusters(int minITSclusters_);
  int getMinITSClusters() const;
  void setMaxITSChi2PerCluster(float maxITSChi2percluster_);
  int getMaxITSChi2PerCluster() const;

  // Calculate logical OR of selection criteria conveniently
  void CombineWithLogicalOR(trackSelectionRequest const& lTraSelRe);

  // Apply track selection checks (to be used in core services)
  template <typename TTrack>
  bool IsTrackSelected(TTrack const& lTrack)
  {
    // Selector that applies all selections
    // Phase-space
    if (lTrack.pt() < minPt)
      return false;
    if (lTrack.pt() > maxPt)
      return false;
    if (lTrack.eta() < minEta)
      return false;
    if (lTrack.eta() > maxEta)
      return false;
    // DCA to PV
    if (std::fabs(lTrack.dcaXY()) < maxDCAz)
      return false;
    // TracksExtra-based
    if (lTrack.hasTPC() == false && requireTPC)
      return false; // FIXME this is a LO approximation
    if (lTrack.tpcNClsFound() < minTPCclusters)
      return false;
    if (lTrack.tpcNClsCrossedRows() < minTPCcrossedrows)
      return false;
    if (lTrack.tpcCrossedRowsOverFindableCls() < minTPCcrossedrowsoverfindable)
      return false;
    if (lTrack.tpcFractionSharedCls() > maxTPCFractionSharedCls)
      return false;
    if (lTrack.hasITS() == false && requireITS)
      return false;
    if (lTrack.itsNCls() < minITSclusters)
      return false;
    if (lTrack.itsChi2NCl() < maxITSChi2percluster)
      return false; // FIXME this is a LO approximation
    return true;
  }
  template <typename TTrack>
  bool IsTrackSelected_TrackExtraCriteria(TTrack const& lTrack)
  {
    // Selector that only applies TracksExtra columns selection
    if (lTrack.hasTPC() == false && requireTPC)
      return false; // FIXME this is a LO approximation
    if (lTrack.tpcNClsFound() < minTPCclusters)
      return false;
    if (lTrack.tpcNClsCrossedRows() < minTPCcrossedrows)
      return false;
    if (lTrack.tpcCrossedRowsOverFindableCls() < minTPCcrossedrowsoverfindable)
      return false;
    if (lTrack.tpcFractionSharedCls() > maxTPCFractionSharedCls)
      return false;
    if (lTrack.hasITS() == false && requireITS)
      return false;
    if (lTrack.itsNCls() < minITSclusters)
      return false;
    if (lTrack.itsChi2NCl() < maxITSChi2percluster)
      return false; // FIXME this is a LO approximation
    return true;
  }

  void SetTightSelections();

  // Helper to print out selections
  void PrintSelections() const;

 private:
  int trackPhysicsType; // 0 - primary, 1 - secondary
  // Phase space (Tracks or TracksIU)
  float minPt;
  float maxPt;
  float minEta;
  float maxEta;
  // DCAs to primary vertex (use for primaries only)
  float maxDCAz;
  float maxDCAxyPtDep;
  // TPC parameters (TracksExtra)
  bool requireTPC; // in Run 3, equiv to hasTPC
  int minTPCclusters;
  int minTPCcrossedrows;
  float minTPCcrossedrowsoverfindable;
  float maxTPCFractionSharedCls;
  // ITS parameters (TracksExtra)
  bool requireITS; // in Run 3, equiv to hasITS
  int minITSclusters;
  float maxITSChi2percluster;

  ClassDefNV(trackSelectionRequest, 3);
};

std::ostream& operator<<(std::ostream& os, trackSelectionRequest const& c);

#endif // COMMON_TOOLS_TRACKSELECTIONREQUEST_H_
