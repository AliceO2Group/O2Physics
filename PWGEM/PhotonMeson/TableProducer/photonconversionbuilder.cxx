// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

/// \file photonconversionbuilder.ch
/// \brief this task produces photon data table with KFParticle.
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>, Tokyo

#include "PWGEM/PhotonMeson/TableProducer/photonconversionbuilder.h"

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  return o2::framework::WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder<o2::aod::V0PhotonsKF, o2::aod::V0Legs, o2::aod::V0LegsXYZ, o2::aod::V0LegsDeDxMC, o2::aod::V0PhotonsPhiVPsi>>(context, o2::framework::TaskName{"photon-conversion-builder"})};
}
