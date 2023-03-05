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

#ifndef TRACKSELECTIONREQUEST_H
#define TRACKSELECTIONREQUEST_H

#include <iosfwd>
#include <Rtypes.h>
#include <TMath.h>

class trackSelectionRequest
{
 public:
  trackSelectionRequest(int minTPCclusters_ = -1, bool option_ = true)
    : minPt{0.0}, maxPt{1e+6}, minEta{-100}, maxEta{+100}, maxDCAz{1e+6}, maxDCAxyPtDep{1e+6}, requireTPCrefit{false}, minTPCclusters{-1}, minTPCcrossedrows{-1}, minTPCcrossedrowsoverfindable{0.0}, requireITSrefit{false}, minITSclusters{-1}, maxITSChi2percluster{1e+6}
  {
    // constructor
  }

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

  void setRequireTPCRefit(bool requireTPCrefit_);
  bool getRequireTPCRefit() const;
  void setMinTPCClusters(int minTPCclusters_);
  int getMinTPCClusters() const;
  void setMinTPCCrossedRows(int minTPCCrossedRows_);
  int getMinTPCCrossedRows() const;
  void setMinTPCCrossedRowsOverFindable(float minTPCCrossedRowsOverFindable_);
  int getMinTPCCrossedRowsOverFindable() const;

  void setRequireITSRefit(bool requireITSrefit_);
  bool getRequireITSRefit() const;
  void setMinITSClusters(int minITSclusters_);
  int getMinITSClusters() const;
  void setMaxITSChi2PerCluster(float maxITSChi2percluster_);
  int getMaxITSChi2PerCluster() const;

  // Calculate logical OR of selection criteria conveniently
  void CombineWithLogicalOR(trackSelectionRequest const& lTraSelRe);

 private:
  // Phase space (Tracks or TracksIU)
  float minPt;
  float maxPt;
  float minEta;
  float maxEta;
  // DCAs to primary vertex (use for primaries only)
  float maxDCAz;
  float maxDCAxyPtDep;
  // TPC parameters (TracksExtra)
  bool requireTPCrefit;
  int minTPCclusters;
  int minTPCcrossedrows;
  float minTPCcrossedrowsoverfindable;
  // ITS parameters (TracksExtra)
  bool requireITSrefit;
  int minITSclusters;
  float maxITSChi2percluster;

  ClassDefNV(trackSelectionRequest, 2);
};

std::ostream& operator<<(std::ostream& os, trackSelectionRequest const& c);

#endif // TRACKSELECTIONREQUEST_H
