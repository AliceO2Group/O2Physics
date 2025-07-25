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
//
// ========================
//
// This code loops over v0 photons and makes pairs for photon HBT analysis.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/PhotonHBT.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/LorentzRotation.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TString.h"

#include <cstring>
#include <iterator>

using namespace o2;
using namespace o2::aod;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonHBT<o2::aod::pwgem::dilepton::core::photonhbt::ggHBTPairType::kPCMEE, MyV0Photons, aod::V0Legs, FilteredMyTracks>>(cfgc, TaskName{"photon-hbt-pcmee"})};
}
