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
// \brief Analyses UD tables (DGCandidates, DGTracks) of DG candidates produced with DGCandProducer
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  01.03.2024

#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/UDGoodRunSelector.h"
#include "PWGUD/Core/decayTree.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct decayTreeAnalyzer {

  // ccdb
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int lastRun = -1;                                          // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB; // bc pattern of colliding bunches

  // goodRun selector
  Configurable<std::string> goodRunsFile{"goodRunsFile", {}, "json with list of good runs"};
  UDGoodRunSelector grsel = UDGoodRunSelector();

  // decay tree object
  Configurable<std::string> parsFile{"parsFile", {}, "json with parameters"};
  decayTree decTree = decayTree();

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // goodRun selector
    grsel.init(goodRunsFile);
    grsel.Print();

    // decay tree object
    decTree.init(parsFile, registry);
    decTree.Print();
  }

  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;

  // PV contributors
  Filter PVContributorFilter = aod::udtrack::isPVContributor == true;
  using PVTracks = soa::Filtered<UDTracksFull>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks, PVTracks const& PVContributors)
  {

    // accept only selected run numbers
    int run = dgcand.runNumber();
    if (!grsel.isGoodRun(run)) {
      return;
    }

    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= 500000) {
      LOGF(info, "Updating bcPattern %d ...", run);
      auto tss = ccdb->getRunDuration(run);
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tss.first);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();
      lastRun = run;
    }

    // is BB bunch?
    auto bcnum = dgcand.globalBC();
    if (run >= 500000 && bcPatternB[bcnum % o2::constants::lhc::LHCMaxBunches] == 0) {
      LOGF(debug, "bcnum[1] %d is not a BB BC", bcnum % o2::constants::lhc::LHCMaxBunches);
      return;
    }

    // check FIT information
    auto bitMin = decTree.dBCRange()[0] + 16;
    auto bitMax = decTree.dBCRange()[1] + 16;
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (decTree.FITvetos()[0] && TESTBIT(dgcand.bbFV0Apf(), bit))
        return;
      if (decTree.FITvetos()[1] && TESTBIT(dgcand.bbFT0Apf(), bit))
        return;
      if (decTree.FITvetos()[2] && TESTBIT(dgcand.bbFT0Cpf(), bit))
        return;
      if (decTree.FITvetos()[3] && TESTBIT(dgcand.bbFDDApf(), bit))
        return;
      if (decTree.FITvetos()[4] && TESTBIT(dgcand.bbFDDCpf(), bit))
        return;
    }

    // check number of PV contributors
    if (dgcand.numContrib() != PVContributors.size()) {
      LOGF(info, "Missmatch of PVContributors %d != %d", dgcand.numContrib(), PVContributors.size());
    }
    auto nTrackRange = decTree.ntrackRange();
    if (dgcand.numContrib() < nTrackRange[0] || dgcand.numContrib() > nTrackRange[1]) {
      LOGF(debug, "Rejected 1: %d not in range [%d, %d].", dgcand.numContrib(), nTrackRange[0], nTrackRange[1]);
      return;
    }

    // skip events with out-of-range rgtrwTOF
    auto rtrwTOF = udhelpers::rPVtrwTOF<true>(dgtracks, PVContributors.size());
    auto minRgtrwTOF = decTree.rgtrTOFMin();
    if (rtrwTOF < minRgtrwTOF) {
      LOGF(debug, "Rejected 3: %f below threshold of %f.", rtrwTOF, minRgtrwTOF);
      return;
    }

    // compute the decay tree
    LOGF(debug, "BC %d", dgcand.globalBC());
    decayTreeResType results = decTree.processTree(PVContributors);
    // decTree.fillHistograms(results, PVContributors);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<decayTreeAnalyzer>(cfgc, TaskName{"decaytree-analyzer"}),
  };
}
