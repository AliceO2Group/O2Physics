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
/// \file Pi0EtaToGammaGammaMCPCMPCMML.cxx
/// \brief This code loops over photons and makes pairs for neutral mesons analyses in MC for PCM-PCM with additional ML-based photon cuts.
/// \author Isabel Kantak, isabel.kantak@cern.ch

#include "PWGEM/PhotonMeson/Core/Pi0EtaToGammaGammaMC.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::pwgem::photonmeson::photonpair;

using MyV0Photons = o2::soa::Filtered<o2::soa::Join<o2::aod::V0PhotonsKF, o2::aod::V0PhotonsPhiVPsi, o2::aod::V0KFEMEventIds, o2::aod::V0PhotonsKFPrefilterBitDerived>>;
using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGammaMC<PairType::kPCMPCM, MyV0Photons, MyMCV0Legs>>(cfgc, TaskName{"pi0eta-to-gammagamma-mc-pcmpcm"}),
  };
}
