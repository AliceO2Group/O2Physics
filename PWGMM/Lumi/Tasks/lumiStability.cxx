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
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsFDD/Digit.h"
#include "Framework/ASoA.h"

using namespace o2;
using namespace o2::framework;
using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

int nBCsPerOrbit = 3564;

namespace o2::aod
{
  namespace trigvertex {
    DECLARE_SOA_COLUMN(TrigVetexCounterFDD, tvcfdd, uint64_t);
    DECLARE_SOA_COLUMN(TrigVetexCounterFT0, tvcft0, uint64_t);
    //DECLARE_SOA_COLUMN(RunNumber, runnumber, uint64_t);
  } //namespace trigvertex

  DECLARE_SOA_TABLE(TrigVertexInfo, "AOD", "TrigVertexInfo",
                    trigvertex::TrigVetexCounterFDD,
                    trigvertex::TrigVetexCounterFT0); //,
                    //trigvertex::RunNumber);
} //namespace o2::aod

struct myExampleTask {
  Produces<o2::aod::TrigVertexInfo> rowtrigvertex;
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisFDDTrigVertex{nBCsPerOrbit + 1, 0.5f, nBCsPerOrbit + 0.5f, "BC in FDD"};
    const AxisSpec axisFT0TrigVertex{nBCsPerOrbit + 1, 0.5f, nBCsPerOrbit + 0.5f, "BC in FT0"};

    // histo about triggers
    histos.add("BCFDD", "BCFDD", kTH1F, {axisFDDTrigVertex});
    histos.add("BCFT0", "BCFT0", kTH1F, {axisFT0TrigVertex});
  }

  void process(aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::BCsWithTimestamps const&)
  {
    int NTrigVertexFDD = 0;
    int NTrigVertexFT0 = 0;

    for( auto const& fdd: fdds){
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if(!bc.timestamp()) continue;

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      std::bitset<8> fddTriggers = fdd.triggerMask();
      bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];

      if(vertex){
        histos.fill(HIST("BCFDD"), localBC);
        NTrigVertexFDD++;
      } // vertex true
    } // loop over FDD events

    for(auto const& ft0 : ft0s){
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if(!bc.timestamp()) continue;

      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      std::bitset<8> fT0Triggers = ft0.triggerMask();
      bool vertex = fT0Triggers[o2::fdd::Triggers::bitVertex];

      if(vertex){
        histos.fill(HIST("BCFT0"), localBC);
        NTrigVertexFT0++;
      } // vertex true
    } // loop over FT0 events
    rowtrigvertex(NTrigVertexFDD,NTrigVertexFT0);
    NTrigVertexFDD = 0;
    NTrigVertexFDD = 0;
  } //end processF
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}