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
// This code is for dielectron analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DileptonHadronMPC.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DileptonHadronMPC<o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron, MyEMH_electron, FilteredMyElectrons>>(cfgc, TaskName{"dielectron-hadron-2pc"})};
}
