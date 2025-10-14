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

/// \file StandardCCDBLoader.cxx
/// \brief A simple object to handle ccdb queries
/// \author ALICE

#ifndef COMMON_TOOLS_STANDARDCCDBLOADER_H_
#define COMMON_TOOLS_STANDARDCCDBLOADER_H_

#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <cmath>
#include <cstdlib>
#include <string>

//__________________________________________
// Standard class to load stuff
// such as matLUT, B and mean Vertex
// partial requests possible.

namespace o2
{
namespace common
{

// ConfigurableGroup with locations
struct StandardCCDBLoaderConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "ccdb";
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
};

class StandardCCDBLoader
{
 public:
  StandardCCDBLoader()
  {
    // constructor - null pointers
    mMeanVtx = nullptr;
    grpmag = nullptr;
    lut = nullptr;
  };

  // commonly needed objects
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  int runNumber = -1;

  template <typename TConfigurableGroup, typename TCCDB, typename TBCs>
  void initCCDBfromBCs(TConfigurableGroup const& cGroup, TCCDB& ccdb, TBCs& bcs, bool getMeanVertex = true)
  {
    // instant load from BCs table. Bonus: protect also against empty bcs
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(cGroup, ccdb, bc.runNumber(), getMeanVertex);
  }

  template <typename TConfigurableGroup, typename TCCDB>
  void initCCDB(TConfigurableGroup const& cGroup, TCCDB& ccdb, int currentRunNumber, bool getMeanVertex = true)
  {
    if (runNumber == currentRunNumber) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << currentRunNumber;
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->template getForRun<o2::base::MatLayerCylSet>(cGroup.lutPath.value, currentRunNumber));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    grpmag = ccdb->template getForRun<o2::parameters::GRPMagField>(cGroup.grpmagPath.value, currentRunNumber);
    LOG(info) << "Setting global propagator magnetic field to current " << grpmag->getL3Current() << " A for run " << currentRunNumber << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    LOG(info) << "Setting global propagator material propagation LUT";
    o2::base::Propagator::Instance()->setMatLUT(lut);
    if (getMeanVertex) {
      // only try this if explicitly requested
      mMeanVtx = ccdb->template getForRun<o2::dataformats::MeanVertexObject>(cGroup.mVtxPath.value, currentRunNumber);
    } else {
      mMeanVtx = nullptr;
    }

    runNumber = currentRunNumber;
  }
};

} // namespace common
} // namespace o2

#endif // COMMON_TOOLS_STANDARDCCDBLOADER_H_
