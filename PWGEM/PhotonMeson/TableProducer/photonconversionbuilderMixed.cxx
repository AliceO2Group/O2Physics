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

/// \file photonconversionbuilderMixed.cxx
/// \brief builder run 2 of the mixed event: rebuilds the photons of the artificial events of photon-mixed-event-producer
///        from their legs (set processRec = false, processMixed = true)
/// \author Stefanie Mrozinski <stefanie.mrozinski@cern.ch>

#include "PWGEM/PhotonMeson/TableProducer/photonconversionbuilder.h"

#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/WorkflowSpec.h>

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  return o2::framework::WorkflowSpec{
    adaptAnalysisTask<PhotonConversionBuilder<o2::aod::V0PhotonsKFMix, o2::aod::V0LegsMix, o2::aod::V0LegsXYZMix, o2::aod::V0LegsDeDxMCMix, o2::aod::V0PhotonsPhiVPsiMix, o2::aod::V0LegsCovMix, o2::aod::V0PhotonsKFChi2Mix, true>>(context, o2::framework::TaskName{"photon-conversion-builder-mixed"})};
}
