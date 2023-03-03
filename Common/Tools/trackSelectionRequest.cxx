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

void trackSelectionRequest::setMinTPCClusters(int minTPCclusters_)
{
  minTPCclusters = minTPCclusters_;
}

int trackSelectionRequest::getMinTPCClusters() const
{
  return minTPCclusters;
}

void trackSelectionRequest::setOption(bool option_)
{
  option = option_;
}

bool trackSelectionRequest::getOption() const
{
  return option;
}
