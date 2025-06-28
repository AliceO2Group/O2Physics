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
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

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

  using HfRedPidDaus2Prong = soa::Join<aod::HfCand2ProngWPid, aod::HfRedPidDau0s_000, aod::HfRedPidDau1s_000>;
  using HfRedPidDaus3Prong = soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfRedPidDau0s_000, aod::HfRedPidDau1s_000, aod::HfRedPidDau2s_000>;

  void process2Prongs(HfRedPidDaus2Prong::iterator const& hfCandPidProngs)
  {
    hfRedPidDau0s(hfCandPidProngs.nSigTpcPi0(), hfCandPidProngs.nSigTofPi0(), hfCandPidProngs.nSigTpcKa0(), hfCandPidProngs.nSigTofKa0(), -999.f, -999.f, hfCandPidProngs.hasTOFProng0(), hfCandPidProngs.hasTPCProng0());
    hfRedPidDau1s(hfCandPidProngs.nSigTpcPi1(), hfCandPidProngs.nSigTofPi1(), hfCandPidProngs.nSigTpcKa1(), hfCandPidProngs.nSigTofKa1(), -999.f, -999.f, hfCandPidProngs.hasTOFProng1(), hfCandPidProngs.hasTPCProng1());
  }
  PROCESS_SWITCH(HfConverterReducedHadronDausPid, process2Prongs, "Produce PID tables for 2-prong candidates", true);

  void process3Prongs(HfRedPidDaus3Prong::iterator const& hfCandPidProngs)
  {
    hfRedPidDau0s(hfCandPidProngs.nSigTpcPi0(), hfCandPidProngs.nSigTofPi0(), hfCandPidProngs.nSigTpcKa0(), hfCandPidProngs.nSigTofKa0(), -999.f, -999.f, hfCandPidProngs.hasTOFProng0(), hfCandPidProngs.hasTPCProng0());
    hfRedPidDau1s(hfCandPidProngs.nSigTpcPi1(), hfCandPidProngs.nSigTofPi1(), hfCandPidProngs.nSigTpcKa1(), hfCandPidProngs.nSigTofKa1(), -999.f, -999.f, hfCandPidProngs.hasTOFProng1(), hfCandPidProngs.hasTPCProng1());
    hfRedPidDau2s(hfCandPidProngs.nSigTpcPi2(), hfCandPidProngs.nSigTofPi2(), hfCandPidProngs.nSigTpcKa2(), hfCandPidProngs.nSigTofKa2(), -999.f, -999.f, hfCandPidProngs.hasTOFProng2(), hfCandPidProngs.hasTPCProng2());
  }
  PROCESS_SWITCH(HfConverterReducedHadronDausPid, process3Prongs, "Produce PID tables for 3-prong candidates", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfConverterReducedHadronDausPid>(cfgc)};
}
