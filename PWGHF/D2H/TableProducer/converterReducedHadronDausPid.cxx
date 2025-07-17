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

/// \file converterReducedHadronDausPid.cxx
/// \brief Task for conversion of daughters pid to version 001
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct HfConverterReducedHadronDausPid {
  Produces<aod::HfRedPidDau0s_001> hfRedPidDau0s;
  Produces<aod::HfRedPidDau1s_001> hfRedPidDau1s;
  Produces<aod::HfRedPidDau2s_001> hfRedPidDau2s;

  using HfRedPidDaus2Prong = soa::Join<aod::HfRed2Prongs, aod::HfRedPidDau0s_000, aod::HfRedPidDau1s_000>;
  using HfRedPidDaus3Prong = soa::Join<aod::HfRed3Prongs, aod::HfRedPidDau0s_000, aod::HfRedPidDau1s_000, aod::HfRedPidDau2s_000>;

  void process2Prongs(HfRedPidDaus2Prong::iterator const& hfCandPidProngs)
  {
    hfRedPidDau0s(hfCandPidProngs.tpcNSigmaPiProng0(), hfCandPidProngs.tofNSigmaPiProng0(), hfCandPidProngs.tpcNSigmaKaProng0(), hfCandPidProngs.tofNSigmaKaProng0(), -999.f, -999.f, hfCandPidProngs.hasTOFProng0(), hfCandPidProngs.hasTPCProng0());
    hfRedPidDau1s(hfCandPidProngs.tpcNSigmaPiProng1(), hfCandPidProngs.tofNSigmaPiProng1(), hfCandPidProngs.tpcNSigmaKaProng1(), hfCandPidProngs.tofNSigmaKaProng1(), -999.f, -999.f, hfCandPidProngs.hasTOFProng1(), hfCandPidProngs.hasTPCProng1());
  }
  PROCESS_SWITCH(HfConverterReducedHadronDausPid, process2Prongs, "Produce PID tables for 2-prong candidates", false);

  void process3Prongs(HfRedPidDaus3Prong::iterator const& hfCandPidProngs)
  {
    hfRedPidDau0s(hfCandPidProngs.tpcNSigmaPiProng0(), hfCandPidProngs.tofNSigmaPiProng0(), hfCandPidProngs.tpcNSigmaKaProng0(), hfCandPidProngs.tofNSigmaKaProng0(), -999.f, -999.f, hfCandPidProngs.hasTOFProng0(), hfCandPidProngs.hasTPCProng0());
    hfRedPidDau1s(hfCandPidProngs.tpcNSigmaPiProng1(), hfCandPidProngs.tofNSigmaPiProng1(), hfCandPidProngs.tpcNSigmaKaProng1(), hfCandPidProngs.tofNSigmaKaProng1(), -999.f, -999.f, hfCandPidProngs.hasTOFProng1(), hfCandPidProngs.hasTPCProng1());
    hfRedPidDau2s(hfCandPidProngs.tpcNSigmaPiProng2(), hfCandPidProngs.tofNSigmaPiProng2(), hfCandPidProngs.tpcNSigmaKaProng2(), hfCandPidProngs.tofNSigmaKaProng2(), -999.f, -999.f, hfCandPidProngs.hasTOFProng2(), hfCandPidProngs.hasTPCProng2());
  }
  PROCESS_SWITCH(HfConverterReducedHadronDausPid, process3Prongs, "Produce PID tables for 3-prong candidates", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfConverterReducedHadronDausPid>(cfgc)};
}
