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

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::string> zvtxsel{"evtflt_zvtx", "zvtx{cwv{rg{-7.0,7.0}-yes:rg{-10.0,10.0}-no,rg{-3.0,3.0}-no}}", "Z vertex cut: zvtx{rg{-7.0,7.0}} or zvtx{cwv{def:var1,var2@,..}}"};
  o2::framework::Configurable<std::string> centmultsel{"evtflt_centmult", "centmult{cwv{mrg{V0M,0,5,10,20,30,40,50,60,70,80}-yes:mrg{CL1,0,5,10,20,30,40,50,60,70,80}-no}}", "Centrality/Multiplicity cut: centmult{mrg{V0M,0,5,10,20,30,40,50,60,70,80}}"};
} eventfilter;

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::string> ttype{"trkflt_ttype", "ttype{FB1-no,FB32-yes,FB64-yes}", "Track types to filter"};
  o2::framework::Configurable<std::string> nclstpc{"trkflt_nclstpc", "nclstpc{cwv{th{70}-yes:th{80}-no,th{90}-no}}", "Min no of TPC clusters: nclstpc{th{70}} or nclstpc{cwv{def:var1,var2,...}}"};
  //  Configurable<std::string> nxrtpc{"trkflt_nxrtpc", "nxrtpc{}", "Min no of TPC crossed rows: nxrtpc{th{70}} or nxrtpc{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> nxrtpc{"trkflt_nxrtpc", "", "Min no of TPC crossed rows: nxrtpc{th{70}} or nxrtpc{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> nclsits{"trkflt_nclsits", "nclsits{}", "Min no of ITS clusters: nclsits{th{3}} or nclsits{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> nclsits{"trkflt_nclsits", "", "Min no of ITS clusters: nclsits{th{3}} or nclsits{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> chi2clustpc{"trkflt_chi2clustpc", "chi2clustpc{cwv{lim{4}-yes:lim{3}-no,lim{90}-no}}", "Max Chi^2 per TPC cluster: chi2clustpc{lim{4}} or chi2clustpc{cwv{def:var1,var2,...}}"};
  //  Configurable<std::string> chi2clusits{"trkflt_chi2clusits", "chi2clusits{}", "Max Chi^2 per ITS cluster: chi2clusits{lim{4}} or chi2clusits{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> chi2clusits{"trkflt_chi2clusits", "", "Max Chi^2 per ITS cluster: chi2clusits{lim{4}} or chi2clusits{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> xrofctpc{"trkflt_xrofctpc", "xrofctpc{}", "Min no of TPC crossed rows over findable clusters: xrofctpc{th{0.70}} or xrofctpc{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> xrofctpc{"trkflt_xrofctpc", "", "Min no of TPC crossed rows over findable clusters: xrofctpc{th{0.70}} or xrofctpc{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> dcaxy{"trkflt_dcaxy", "dcaxy{}", "Max DCAxy: dcaxy{lim{2.3}} or dcaxy{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> dcaxy{"trkflt_dcaxy", "", "Max DCAxy: dcaxy{lim{2.3}} or dcaxy{cwv{def;var1,var2,...}}"};
  //  Configurable<std::string> dcaz{"trkflt_dcaz", "dcaz{}", "Max DCAz: dcaz{lim{3.0}} or dcaz{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> dcaz{"trkflt_dcaz", "", "Max DCAz: dcaz{lim{3.0}} or dcaz{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> ptrange{"trkflt_pt", "pT{rg{0.2,10.0}}", "pT range: pT{th{0.2}} or pT{cwv{def;var1,var2,...}}"};
  o2::framework::Configurable<std::string> etarange{"trkflt_eta", "eta{rg{-0.8,0.8}}", "eta range: eta{rg{-0.9,0.9}} or eta{cwv{def;var1,var2,...}}"};
} trackfilter;

struct : o2::framework::ConfigurableGroup {
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::string> itsel{"pidflt_its_el", "itsel{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the its electron dE/dx line, below/above"};
    o2::framework::Configurable<std::string> itsmu{"pidflt_its_mu", "itsmu{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the its muon dE/dx line, below/above"};
    o2::framework::Configurable<std::string> itspi{"pidflt_its_pi", "itspi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the its pion dE/dx line, below/above"};
    o2::framework::Configurable<std::string> itska{"pidflt_its_ka", "itska{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the its kaon dE/dx line, below/above"};
    o2::framework::Configurable<std::string> itspr{"pidflt_its_pr", "itspr{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the its proton dE/dx line, below/above"};
  } piditsfilter;
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::string> tpcel{"pidflt_tpc_el", "tpcel{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tpc electron dE/dx line, below/above"};
    o2::framework::Configurable<std::string> tpcmu{"pidflt_tpc_mu", "tpcmu{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tpc muon dE/dx line, below/above"};
    o2::framework::Configurable<std::string> tpcpi{"pidflt_tpc_pi", "tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tpc pion dE/dx line, below/above"};
    o2::framework::Configurable<std::string> tpcka{"pidflt_tpc_ka", "tpcks{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tpc kaon dE/dx line, below/above"};
    o2::framework::Configurable<std::string> tpcpr{"pidflt_tpc_pr", "tpcpr{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tpc proton dE/dx line, below/above"};
  } pidtpcfilter;
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::string> tpcel{"pidflt_tof_el", "tofel{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tof electron line, below/above"};
    o2::framework::Configurable<std::string> tpcmu{"pidflt_tof_mu", "tofmu{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tof muon line, below/above"};
    o2::framework::Configurable<std::string> tpcpi{"pidflt_tof_pi", "tofpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tof pion line, below/above"};
    o2::framework::Configurable<std::string> tpcka{"pidflt_tof_ka", "tofks{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tof kaon line, below/above"};
    o2::framework::Configurable<std::string> tpcpr{"pidflt_tof_pr", "tofpr{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}", "nsigmas to the tof proton line, below/above"};
  } pidtoffilter;
} pidfilter;

#endif // O2_ANALYSIS_DPTDPTSKIMCONF_H
