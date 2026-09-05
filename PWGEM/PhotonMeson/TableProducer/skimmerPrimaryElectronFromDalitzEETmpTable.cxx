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

/// \file skimmerPrimaryElectronFromDalitzEETmpTable.cxx
/// \brief write electrons for Dalitz into a temporary table that is then used for event skimming
/// \author joshua.konig@cern.ch

#include "PWGEM/PhotonMeson/TableProducer/skimmerPrimaryElectronFromDalitzEE.h"

using namespace o2::framework;
using namespace o2::aod;
using namespace o2::framework::expressions;

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(context);

  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectronFromDalitzEE<o2::aod::EMPrimaryElectronsFromDalitzTmp, o2::aod::EMPrimaryElectronsDeDxMCTmp, o2::aod::V0PhotonsKFTmp>>(context, TaskName{"skimmer-primary-electron-from-dalitzee-tmptable"})};
}
