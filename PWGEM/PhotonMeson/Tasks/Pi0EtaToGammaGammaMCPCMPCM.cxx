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
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/Pi0EtaToGammaGammaMC.h"

using namespace o2;
using namespace o2::aod;

// using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
// using MyV0Photon = MyV0Photons::iterator;

// using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
// using MyMCV0Leg = MyMCV0Legs::iterator;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGammaMC<PairType::kPCMPCM, MyV0Photons, MyMCV0Legs>>(cfgc, TaskName{"pi0eta-to-gammagamma-mc-pcmpcm"}),
  };
}
