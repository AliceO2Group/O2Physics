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
// This code loops over photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGEM/PhotonMeson/Core/TaggingPi0.h"

using namespace o2;
using namespace o2::aod;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TaggingPi0<PairType::kPCMDalitzEE, MyV0Photons, aod::V0Legs, MyPrimaryElectrons>>(cfgc, TaskName{"tagging-pi0-pcmdalitzee"}),
  };
}
