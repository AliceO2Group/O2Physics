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
// Analysis task for single electron QC in MC
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "PWGEM/Dilepton/Core/SingleTrackQCMC.h"

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SingleTrackQCMC<o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron, FilteredMyMCElectrons, MySmearedElectrons>>(cfgc, TaskName{"single-electron-qc-mc"})};
}
