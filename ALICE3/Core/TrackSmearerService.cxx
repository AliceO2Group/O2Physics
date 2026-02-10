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

/// \file TrackSmearerService.cxx
/// \brief Implementation for smearer service for the on-the-fly simulation
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>


#include "ALICE3/Core/TrackSmearerService.h"
#include <Framework/CommonServices.h>
#include <Framework/Plugins.h>
#include <Framework/ServiceHandle.h>
#include <Framework/ServiceSpec.h>

#include <string>
#include <TPDGCode.h>

#include <vector>

bool o2::upgrade::TrackSmearerImpl::mIsInit = false;
bool o2::upgrade::TrackSmearerImpl::mCleanLutWhenLoaded = true;

void o2::upgrade::TrackSmearerImpl::initSmearer(o2::ccdb::BasicCCDBManager* ccdb, bool cleanLutWhenLoaded)
{
  if (mIsInit) {
    LOG(fatal) << "TrackSmearerImpl already initialized, cannot re-initialize";
  }

  if (!ccdb) {
    LOG(fatal) << "CCDB manager is not set, cannot initialize TrackSmearerImpl";
  }

  mIsInit = true;
  mCcdb = ccdb;
  mCcdb->setURL("http://alice-ccdb.cern.ch");
  mCcdb->setTimestamp(-1);
  mCleanLutWhenLoaded = cleanLutWhenLoaded;
}


std::vector<std::unique_ptr<o2::delphes::DelphesO2TrackSmearer>> o2::upgrade::TrackSmearerImpl::smearerContainer;
void o2::upgrade::TrackSmearerImpl::initCfg(int icfg, std::map<std::string, std::string> globalConfiguration)
{
  smearerContainer.emplace_back(std::make_unique<o2::delphes::DelphesO2TrackSmearer>());
  smearerContainer[icfg]->setCcdbManager(mCcdb);
  smearerContainer[icfg]->setDownloadPath("./.ALICE3/");
  smearerContainer[icfg]->setCleanupDownloadedFile(mCleanLutWhenLoaded);
  for (const auto& entry : globalConfiguration) {
    int pdg = 0;
    if (entry.first.find("lut") != 0) {
      continue;
    }
    if (entry.first.find("lutEl") != std::string::npos) {
      pdg = kElectron;
    } else if (entry.first.find("lutMu") != std::string::npos) {
      pdg = kMuonMinus;
    } else if (entry.first.find("lutPi") != std::string::npos) {
      pdg = kPiPlus;
    } else if (entry.first.find("lutKa") != std::string::npos) {
      pdg = kKPlus;
    } else if (entry.first.find("lutPr") != std::string::npos) {
      pdg = kProton;
    } else if (entry.first.find("lutDe") != std::string::npos) {
      pdg = o2::constants::physics::kDeuteron;
    } else if (entry.first.find("lutTr") != std::string::npos) {
      pdg = o2::constants::physics::kTriton;
    } else if (entry.first.find("lutHe3") != std::string::npos) {
      pdg = o2::constants::physics::kHelium3;
    } else if (entry.first.find("lutAl") != std::string::npos) {
      pdg = o2::constants::physics::kAlpha;
    }

    std::string filename = entry.second;
    if (pdg == 0) {
      LOG(fatal) << "Unknown LUT entry " << entry.first << " for global configuration";
    }
    LOG(info) << "Loading LUT for pdg " << pdg << " for config " << icfg << " from provided file '" << filename << "'";
    if (filename.empty()) {
      LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
    }
    // strip from leading/trailing spaces
    filename.erase(0, filename.find_first_not_of(" "));
    filename.erase(filename.find_last_not_of(" ") + 1);
    if (filename.empty()) {
      LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
    }
    bool success = smearerContainer[icfg]->loadTable(pdg, filename.c_str());
    if (!success) {
      LOG(fatal) << "Having issue with loading the LUT " << pdg << " " << filename;
    }
    LOG(info) << "Successfully loaded LUT for pdg " << pdg;
  }
}

o2::delphes::DelphesO2TrackSmearer* o2::upgrade::TrackSmearerImpl::getSmearer(int icfg)
{
  if (icfg < 0 || icfg >= static_cast<int>(smearerContainer.size())) {
    LOG(fatal) << "Configuration index " << icfg << " out of bounds";
    return nullptr;
  }
  if (!smearerContainer[icfg]) {
    LOG(fatal) << "nullptr entry in tracksmearer container";
    return nullptr;
  }
  return smearerContainer[icfg].get();
}


void o2::upgrade::TrackSmearerImpl::setReady()
{
  std::ofstream okFile(".TrackSmearerOK");
  okFile.close();
  LOG(info) << "Track smearer ready";
}

void o2::upgrade::TrackSmearerImpl::waitReady()
{
  while (!std::ifstream(".TrackSmearerOK")) {
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
}

struct TrackSmearerSupport : o2::framework::ServicePlugin {
  o2::framework::ServiceSpec* create() final
  {
    return new o2::framework::ServiceSpec{
      .name = "track-smearer",
      .init = [](o2::framework::ServiceRegistryRef, o2::framework::DeviceState&, fair::mq::ProgOptions&) -> o2::framework::ServiceHandle {
        auto* wrapper = new o2::upgrade::TrackSmearerContainer();
        auto* ptr = new o2::upgrade::TrackSmearerImpl();
        wrapper->setInstance(ptr);
        return o2::framework::ServiceHandle{o2::framework::TypeIdHelpers::uniqueId<o2::upgrade::TrackSmearerContainer>(),
                                          wrapper,
                                          o2::framework::ServiceKind::Serial,
                                          "database-pdg"};
      },
      .configure = o2::framework::CommonServices::noConfiguration(),
      .exit = [](o2::framework::ServiceRegistryRef, void* service) {
        auto* s = reinterpret_cast<o2::upgrade::TrackSmearerContainer*>(service);
        delete s;
      },
      .kind = o2::framework::ServiceKind::Serial
    };
  }
};

DEFINE_DPL_PLUGINS_BEGIN
DEFINE_DPL_PLUGIN_INSTANCE(TrackSmearerSupport, CustomService);
DEFINE_DPL_PLUGINS_END

