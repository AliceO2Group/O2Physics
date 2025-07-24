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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GENERICFRAMEWORKLINKDEF_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GENERICFRAMEWORKLINKDEF_H_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class GFWPowerArray + ;
#pragma link C++ class GFWCumulant + ;
#pragma link C++ class GFW + ;
#pragma link C++ class ProfileSubset + ;
#pragma link C++ class FlowContainer + ;
#pragma link C++ class GFWWeights + ;
#pragma link C++ class GFWWeightsList + ;
#pragma link C++ class BootstrapProfile + ;
#pragma link C++ class FlowPtContainer + ;
#pragma link C++ class o2::analysis::genericframework::GFWBinningCuts + ;
#pragma link C++ class o2::analysis::genericframework::GFWRegions + ;
#pragma link C++ class o2::analysis::genericframework::GFWCorrConfigs + ;

#endif // PWGCF_GENERICFRAMEWORK_CORE_GENERICFRAMEWORKLINKDEF_H_
