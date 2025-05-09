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

#ifndef COMMON_TOOLS_STANDARDCCDBLOADER_H_
#define COMMON_TOOLS_STANDARDCCDBLOADER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "Framework/AnalysisDataModel.h"

//__________________________________________
// Standard class to load stuff
// such as matLUT, B and mean Vertex
// partial requests possible.

class StandardCCDBLoader
{
 public:
  StandardCCDBLoader()
  {
    // constructor - sets standard locations
    ccdburl = "http://alice-ccdb.cern.ch";
    lutPath = "GLO/Param/MatLUT";
    geoPath = "GLO/Config/GeometryAligned";
    grpmagPath = "GLO/Config/GRPMagField";
    mVtxPath = "GLO/Calib/MeanVertex";

    mMeanVtx = nullptr;
    grpmag = nullptr;
    lut = nullptr;
  };

  // commonly needed objects
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  std::string ccdburl;
  std::string lutPath;
  std::string geoPath;
  std::string grpmagPath;
  std::string mVtxPath;
  int runNumber = -1;

  // takes a configurableGroup and reads in necessary configs
  template <typename TConfigurableGroup>
  void readConfiguration(TConfigurableGroup const& cGroup)
  {
    ccdburl = cGroup.ccdburl.value;
    lutPath = cGroup.lutPath.value;
    geoPath = cGroup.geoPath.value;
    grpmagPath = cGroup.grpmagPath.value;
    mVtxPath = cGroup.mVtxPath.value;
    LOGF(info, "Standard CCDB loader configuration acquired.");
  }

  template <typename TCCDB, typename TBCs>
  void initCCDBfromBCs(TCCDB& ccdb, TBCs& bcs, bool getMeanVertex = true)
  {
    // instant load from BCs table. Bonus: protect also against empty bcs
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(ccdb, bc.runNumber(), getMeanVertex);
  }

  template <typename TCCDB>
  void initCCDB(TCCDB& ccdb, int currentRunNumber, bool getMeanVertex = true)
  {
    if (runNumber == currentRunNumber) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << currentRunNumber;
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->template getForRun<o2::base::MatLayerCylSet>(lutPath, currentRunNumber));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    grpmag = ccdb->template getForRun<o2::parameters::GRPMagField>(grpmagPath, currentRunNumber);
    LOG(info) << "Setting global propagator magnetic field to current " << grpmag->getL3Current() << " A for run " << currentRunNumber << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    LOG(info) << "Setting global propagator material propagation LUT";
    o2::base::Propagator::Instance()->setMatLUT(lut);
    if (getMeanVertex) {
      // only try this if explicitly requested
      mMeanVtx = ccdb->template getForRun<o2::dataformats::MeanVertexObject>(mVtxPath, currentRunNumber);
    } else {
      mMeanVtx = nullptr;
    }

    runNumber = currentRunNumber;
  }
};

#endif // COMMON_TOOLS_STANDARDCCDBLOADER_H_
