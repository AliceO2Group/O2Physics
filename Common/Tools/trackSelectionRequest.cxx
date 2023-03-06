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

// see header for a more detailed description.

#include "trackSelectionRequest.h"
#include <iostream>

std::ostream& operator<<(std::ostream& os, trackSelectionRequest const& c)
{
  os << "MinTPCClusters value: " << c.getMinTPCClusters();
  return os;
}

void trackSelectionRequest::setMinPt(float minPt_)
{
  minPt = minPt_;
}
int trackSelectionRequest::getMinPt() const
{
  return minPt;
}
void trackSelectionRequest::setMaxPt(float maxPt_)
{
  maxPt = maxPt_;
}
int trackSelectionRequest::getMaxPt() const
{
  return maxPt;
}
void trackSelectionRequest::setMinEta(float minEta_)
{
  minEta = minEta_;
}
int trackSelectionRequest::getMinEta() const
{
  return minEta;
}
void trackSelectionRequest::setMaxEta(float maxEta_)
{
  maxEta = maxEta_;
}
int trackSelectionRequest::getMaxEta() const
{
  return maxEta;
}
void trackSelectionRequest::setMaxDCAz(float maxDCAz_)
{
  maxDCAz = maxDCAz_;
}
int trackSelectionRequest::getMaxDCAz() const
{
  return maxDCAz;
}
void trackSelectionRequest::setMaxDCAxyPtDep(float maxDCAxyPtDep_)
{
  maxDCAxyPtDep = maxDCAxyPtDep_;
}
int trackSelectionRequest::getMaxDCAxyPtDep() const
{
  return maxDCAxyPtDep;
}
void trackSelectionRequest::setRequireTPCRefit(bool requireTPCrefit_)
{
  requireTPCrefit = requireTPCrefit_;
}
bool trackSelectionRequest::getRequireTPCRefit() const
{
  return requireTPCrefit;
}
void trackSelectionRequest::setMinTPCClusters(int minTPCclusters_)
{
  minTPCclusters = minTPCclusters_;
}
int trackSelectionRequest::getMinTPCClusters() const
{
  return minTPCclusters;
}
void trackSelectionRequest::setMinTPCCrossedRows(int minTPCcrossedrows_)
{
  minTPCcrossedrows = minTPCcrossedrows_;
}
int trackSelectionRequest::getMinTPCCrossedRows() const
{
  return minTPCcrossedrows;
}
void trackSelectionRequest::setMinTPCCrossedRowsOverFindable(float minTPCcrossedrowsoverfindable_)
{
  minTPCcrossedrowsoverfindable = minTPCcrossedrowsoverfindable_;
}
int trackSelectionRequest::getMinTPCCrossedRowsOverFindable() const
{
  return minTPCcrossedrowsoverfindable;
}
void trackSelectionRequest::setRequireITSRefit(bool requireITSrefit_)
{
  requireITSrefit = requireITSrefit_;
}
bool trackSelectionRequest::getRequireITSRefit() const
{
  return requireITSrefit;
}
void trackSelectionRequest::setMinITSClusters(int minITSclusters_)
{
  minITSclusters = minITSclusters_;
}
int trackSelectionRequest::getMinITSClusters() const
{
  return minITSclusters;
}
void trackSelectionRequest::setMaxITSChi2PerCluster(float maxITSChi2percluster_)
{
  maxITSChi2percluster = maxITSChi2percluster_;
}
int trackSelectionRequest::getMaxITSChi2PerCluster() const
{
  return maxITSChi2percluster;
}

void trackSelectionRequest::CombineWithLogicalOR(trackSelectionRequest const& lTraSelRe)
{
  // This helper method provides the ability to conveniently combine
  // several sets of track selection requests
  if (lTraSelRe.getMinPt() < minPt)
    minPt = lTraSelRe.getMinPt();
  if (lTraSelRe.getMaxPt() > maxPt)
    maxPt = lTraSelRe.getMaxPt();
  if (lTraSelRe.getMinEta() < minEta)
    minEta = lTraSelRe.getMinEta();
  if (lTraSelRe.getMaxEta() > maxEta)
    maxEta = lTraSelRe.getMaxEta();

  if (lTraSelRe.getMaxDCAz() > maxDCAz)
    maxDCAz = lTraSelRe.getMaxDCAz();
  if (lTraSelRe.getMaxDCAxyPtDep() > maxDCAxyPtDep)
    maxDCAxyPtDep = lTraSelRe.getMaxDCAxyPtDep();

  if (lTraSelRe.getRequireTPCRefit() == false)
    requireTPCrefit = false;
  if (lTraSelRe.getMinTPCClusters() < minTPCclusters)
    minTPCclusters = lTraSelRe.getMinTPCClusters();
  if (lTraSelRe.getMinTPCCrossedRows() < minTPCcrossedrows)
    minTPCcrossedrows = lTraSelRe.getMinTPCCrossedRows();
  if (lTraSelRe.getMinTPCCrossedRowsOverFindable() < minTPCcrossedrowsoverfindable)
    minTPCcrossedrowsoverfindable = lTraSelRe.getMinTPCCrossedRowsOverFindable();

  if (lTraSelRe.getRequireITSRefit() == false)
    requireITSrefit = false;
  if (lTraSelRe.getMinITSClusters() < minITSclusters)
    minITSclusters = lTraSelRe.getMinITSClusters();
  if (lTraSelRe.getMaxITSChi2PerCluster() > maxITSChi2percluster)
    maxITSChi2percluster = lTraSelRe.getMaxITSChi2PerCluster();
  return;
}