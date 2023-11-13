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

#ifndef PWGMM_HISTOGRAMS_H
#define PWGMM_HISTOGRAMS_H
#include "Axes.h"

namespace pwgmm::mult {
using namespace o2::framework;
static constexpr std::array<std::array<std::string_view, 2>, 2> categories{
  {
   {"Tracks", "Events"},                       //
   {"Tracks/Centrality", "Events/Centrality"}  //
  }                                            //
};
}

#endif // PWGMM_HISTOGRAMS_H
