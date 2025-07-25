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

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"

#include <string>
#include <vector>

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<std::string>> bfield{"evtflt_bfield", {"positive-yes", "negative-yes"}, "B filed polarity cut: both 'yes' default, anything else alternative"};
  o2::framework::Configurable<std::vector<std::string>> zvtxsel{"evtflt_zvtx", {"rg{-7.0,7.0}-yes", "rg{-10.0,10.0}-no", "rg{-3.0,3.0}-no"}, "Z vertex cut: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> centmultsel{"evtflt_centmult", {"mrg{V0M,0,5,10,20,30,40,50,60,70,80}-yes", "mrg{CL1,0,5,10,20,30,40,50,60,70,80}-no"}, "Centrality/Multiplicity cut: first, default, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> pileuprej{"evtflt_pileuprej", {"fnrg{V0M_TPCOUT=-0.5+3.7*x-0.14*x*x,0.5-3.7*x+0.14*x*x}-yes", "fnrg{V0M_TRKLETS=-0.5+3.7*x-0.14*x*x,0.5-3.7*x+0.14*x*x}-no"}, "Advanced pile-up rejection cut: first, default, next, alternatives"};
} eventfilter;

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<std::string>> ttype{"trkflt_ttype", {"FB1-no", "FB32-yes", "FB64-yes"}, "Track types to filter"};
  o2::framework::Configurable<std::vector<std::string>> nclstpc{"trkflt_nclstpc", {"th{70}-yes", "th{80}-no", "th{90}-no"}, "Min no of TPC clusters: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> nxrtpc{"trkflt_nxrtpc", {""}, "Min no of TPC crossed rows: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> nclsits{"trkflt_nclsits", {""}, "Min no of ITS clusters: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> chi2clustpc{"trkflt_chi2clustpc", {"lim{4}-yes", "lim{3}-no", "lim{90}-no"}, "Max Chi^2 per TPC cluster: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> chi2clusits{"trkflt_chi2clusits", {""}, "Max Chi^2 per ITS cluster: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> xrofctpc{"trkflt_xrofctpc", {""}, "Min no of TPC crossed rows over findable clusters: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> dcaxy{"trkflt_dcaxy", {""}, "Max DCAxy: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> dcaz{"trkflt_dcaz", {""}, "Max DCAz: first, default value, next, alternatives"};
  o2::framework::Configurable<std::vector<std::string>> ptrange{"trkflt_pt", {"rg{0.1,50.0}"}, "pT range"};
  o2::framework::Configurable<std::vector<std::string>> etarange{"trkflt_eta", {"rg{-1.0,1.0}"}, "eta range"};
} trackfilter;

struct : o2::framework::ConfigurableGroup {
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::vector<std::string>> tpcel{"pidflt_tpc_el", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tpc electron dE/dx line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcmu{"pidflt_tpc_mu", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tpc muon dE/dx line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcpi{"pidflt_tpc_pi", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tpc pion dE/dx line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcka{"pidflt_tpc_ka", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tpc kaon dE/dx line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcpr{"pidflt_tpc_pr", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tpc proton dE/dx line, below/above: first, default value, next, alternatives"};
  } pidtpcfilter;
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::vector<std::string>> tpcel{"pidflt_tof_el", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tof electron line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcmu{"pidflt_tof_mu", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tof muon line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcpi{"pidflt_tof_pi", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tof pion line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcka{"pidflt_tof_ka", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tof kaon line, below/above: first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> tpcpr{"pidflt_tof_pr", {"rg{-3.0,3.0}-yes", "rg{-2.0,2.0}-no", "rg{-3.0,5.0}-no"}, "nsigmas to the tof proton line, below/above: first, default value, next, alternatives"};
  } pidtoffilter;
#ifdef INCORPORATEBAYESIANPID
  struct : ConfigurableGroup {
    o2::framework::Configurable<std::vector<std::string>> bayel{"pidflt_bayes_el", {"th{80}-yes", "th{70}-no", "th{90}-no"}, "Bayesian probability for electron (%%): first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> baymu{"pidflt_bayes_mu", {"th{80}-yes", "th{70}-no", "th{90}-no"}, "Bayesian probability for muon (%%): first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> baypi{"pidflt_bayes_pi", {"th{80}-yes", "th{70}-no", "th{90}-no"}, "Bayesian probability for pion (%%): first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> bayka{"pidflt_bayes_ka", {"th{80}-yes", "th{70}-no", "th{90}-no"}, "Bayesian probability for kaon (%%): first, default value, next, alternatives"};
    o2::framework::Configurable<std::vector<std::string>> baypr{"pidflt_bayes_pr", {"th{80}-yes", "th{70}-no", "th{90}-no"}, "Bayesian probability for proton (%%): first, default value, next, alternatives"};
  } pidbayesfilter;
#endif
} pidfilter;

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://ccdb-test.cern.ch:8080", "url of the skimming ccdb repository"};
  o2::framework::Configurable<std::string> ccdbpath{"ccdbpath", "Users/v/victor/Skimming", "url of the skimming ccdb repository"};
  o2::framework::Configurable<std::string> filterdate{"filterdate", "20221115", "the date for the skimming production with the current filter configuration"};
} filterccdb;
