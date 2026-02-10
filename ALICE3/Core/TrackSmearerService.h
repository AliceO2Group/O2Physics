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

/// \file TrackSmearerService.h
/// \brief Implementation for smearer service for the on-the-fly simulation
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>


#ifndef ALICE3_CORE_TRACKSMEARERSERVICE_H_
#define ALICE3_CORE_TRACKSMEARERSERVICE_H_

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "Framework/Plugins.h"
#include "CCDB/BasicCCDBManager.h"
#include "ALICE3/Core/FastTracker.h"

#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>

namespace o2::upgrade {

struct TrackSmearerImpl {
  static std::vector<std::unique_ptr<o2::delphes::DelphesO2TrackSmearer>> smearerContainer;
  void initSmearer(o2::ccdb::BasicCCDBManager* ccdb, bool cleanLutWhenLoaded);
  void initCfg(int icfg, std::map<std::string, std::string> globalConfiguration);
  size_t size() { return smearerContainer.size(); }
  void setReady();
  void waitReady();
  o2::delphes::DelphesO2TrackSmearer* getSmearer(int icfg = 0);


private:
  static bool mIsInit;
  static bool mCleanLutWhenLoaded;
  o2::ccdb::BasicCCDBManager* mCcdb = nullptr;
};

struct TrackSmearerContainer : o2::framework::LoadableServicePlugin<TrackSmearerImpl> {
  TrackSmearerContainer() : LoadableServicePlugin{"O2PhysicsALICE3Core:TrackSmearerSupport"}
  {
  }
};

} // namespace o2::upgrade

#endif // ALICE3_CORE_TRACKSMEARERSERVICE_H_