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

static constexpr double default_matrix[3][3] = {{1.1, 1.2, 1.3}, {2.1, 2.2, 2.3}, {3.1, 3.2, 3.3}};

class trackSelectionRequest
{
 public:
  trackSelectionRequest(int minTPCclusters_ = 70, bool option_ = true)
    : minTPCclusters{minTPCclusters_}, option{option_}
  {
  }

  void setMinTPCClusters(int minTPCclusters_);
  int getMinTPCClusters() const;

  void setOption(bool option_);
  bool getOption() const;

 private:
  int minTPCclusters;
  bool option;

  ClassDefNV(trackSelectionRequest, 1);
};

std::ostream& operator<<(std::ostream& os, trackSelectionRequest const& c);

#endif // TRACKSELECTIONREQUEST_H
