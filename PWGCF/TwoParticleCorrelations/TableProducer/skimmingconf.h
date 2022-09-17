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
#ifndef O2_ANALYSIS_DPTDPTSKIMCONF_H
#define O2_ANALYSIS_DPTDPTSKIMCONF_H

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

struct : ConfigurableGroup {
  Configurable<std::string> zvtxsel{"evtflt_zvtx", "zvtx{cwv{rg{-7.0,7.0}:rg{-10.0,10.0}-no,rg{-3.0,3.0}-no}}", "Z vertex cut: zvtx{rg{-7.0,7.0}} or zvtx{cwv{def:var1,var2@,..}}"};
  Configurable<std::string> centmultsel{"evtflt_centmult", "centmult{cwv{mrg{V0M,0,5,10,20,30,40,50,60,70,80}-yes:mrg{CL1,0,5,10,20,30,40,50,60,70,80}-no}}", "Centrality/Multiplicity cut: centmult{mrg{V0M,0,5,10,20,30,40,50,60,70,80}}"};
} eventfilter;

struct : ConfigurableGroup {
  Configurable<std::string> ttype{"trkflt_ttype", "ttype{FB1,FB32,FB64}", "Track types to filter"};
  Configurable<std::string> nclstpc{"trkflt_nclstpc", "nclstpc{cwv{th{70}:th{80}-no,th{90}-no}}", "Min no of TPC clusters: nclstpc{th{70}} or nclstpc{cwv{def:var1,var2,...}}"};
  //  Configurable<std::string> nxrtpc{"trkflt_nxrtpc", "nxrtpc{}", "Min no of TPC crossed rows: nxrtpc{th{70}} or nxrtpc{cwv{def;var1,var2,...}}"};
  Configurable<std::string> nxrtpc{"trkflt_nxrtpc", "", "Min no of TPC crossed rows: nxrtpc{th{70}} or nxrtpc{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> nclsits{"trkflt_nclsits", "nclsits{}", "Min no of ITS clusters: nclsits{th{3}} or nclsits{cwv{def;var1,var2,...}}"};
  Configurable<std::string> nclsits{"trkflt_nclsits", "", "Min no of ITS clusters: nclsits{th{3}} or nclsits{cwv{def;var1,var2,...}}"};
  Configurable<std::string> chi2clustpc{"trkflt_chi2clustpc", "chi2clustpc{cwv{lim{4}:lim{3}-no,lim{90}-no}}", "Max Chi^2 per TPC cluster: chi2clustpc{lim{4}} or chi2clustpc{cwv{def:var1,var2,...}}"};
  //  Configurable<std::string> chi2clusits{"trkflt_chi2clusits", "chi2clusits{}", "Max Chi^2 per ITS cluster: chi2clusits{lim{4}} or chi2clusits{cwv{def;var1,var2,...}}"};
  Configurable<std::string> chi2clusits{"trkflt_chi2clusits", "", "Max Chi^2 per ITS cluster: chi2clusits{lim{4}} or chi2clusits{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> xrofctpc{"trkflt_xrofctpc", "xrofctpc{}", "Min no of TPC crossed rows over findable clusters: xrofctpc{th{0.70}} or xrofctpc{cwv{def;var1,var2,...}}"};
  Configurable<std::string> xrofctpc{"trkflt_xrofctpc", "", "Min no of TPC crossed rows over findable clusters: xrofctpc{th{0.70}} or xrofctpc{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> dcaxy{"trkflt_dcaxy", "dcaxy{}", "Max DCAxy: dcaxy{lim{2.3}} or dcaxy{cwv{def;var1,var2,...}}"};
  Configurable<std::string> dcaxy{"trkflt_dcaxy", "", "Max DCAxy: dcaxy{lim{2.3}} or dcaxy{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> dcaz{"trkflt_dcaz", "dcaz{}", "Max DCAz: dcaz{lim{3.0}} or dcaz{cwv{def;var1,var2,...}}"};
  Configurable<std::string> dcaz{"trkflt_dcaz", "", "Max DCAz: dcaz{lim{3.0}} or dcaz{cwv{def;var1,var2,...}}"};
  Configurable<std::string> ptrange{"trkflt_pt", "pT{rg{0.2,10.0}}", "pT range: pT{th{0.2}} or pT{cwv{def;var1,var2,...}}"};
  Configurable<std::string> etarange{"trkflt_eta", "eta{rg{-0.8,0.8}}", "eta range: eta{rg{-0.9,0.9}} or eta{cwv{def;var1,var2,...}}"};
} trackfilter;

struct : ConfigurableGroup {

} pidfilter;

#endif // O2_ANALYSIS_DPTDPTSKIMCONF_H
