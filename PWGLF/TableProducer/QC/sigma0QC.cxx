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
// This is a task that employs the standard derived V0 tables and attempts to combine
// two V0s into a Sigma0 -> Lambda + gamma candidate. 
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Sigma0 QC task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "CCDB/BasicCCDBManager.h"
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <PWGLF/Utils/sigma0BuilderHelper.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using StraColls = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>;


struct sigma0Sorted {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  // Load configurables and the sigma0 module, please
  o2::pwglf::sigma0::evselConfigurables evselOpts;
  o2::pwglf::sigma0::lambdaselConfigurables lambdaselOpts;
  o2::pwglf::sigma0::photonselConfigurables photonselOpts;
  o2::pwglf::sigma0::sigma0selConfigurables sigma0selOpts;
  o2::pwglf::sigma0::pi0selConfigurables pi0selOpts;
  o2::pwglf::sigma0::axisConfigurables axisOpts;

  o2::pwglf::sigma0::Sigma0BuilderModule sigma0Module;

  // For manual sliceBy  
  SliceCache cache;
  Preslice<V0StandardDerivedDatas> perCollisionSTDDerived = o2::aod::v0data::straCollisionId;  

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> fverbose{"fverbose", true, "QA printout."};

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Initialize task
    sigma0Module.init(histos, evselOpts, lambdaselOpts, photonselOpts, sigma0selOpts, pi0selOpts, axisOpts);
  }

  // Dummy process function 
  void process(StraColls const&){} 

  void processRealData(StraColls const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&){
   auto start = std::chrono::high_resolution_clock::now();

   sigma0Module.process(collisions, fullV0s, histos, cache, ccdb, rateFetcher); 
   
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;

   if (fverbose) LOGF(info, "[Process function call, Sorted] N. Collisions: %i, N. V0s: %i, Processing time (s): %lf", collisions.size(), fullV0s.size(), elapsed.count());      
  }
  
  PROCESS_SWITCH(sigma0Sorted, processRealData, "process run 3 real data", true);    
};

struct sigma0Unsorted {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  // Load configurables and the sigma0 module, please
  o2::pwglf::sigma0::evselConfigurables evselOpts;
  o2::pwglf::sigma0::lambdaselConfigurables lambdaselOpts;
  o2::pwglf::sigma0::photonselConfigurables photonselOpts;
  o2::pwglf::sigma0::sigma0selConfigurables sigma0selOpts;
  o2::pwglf::sigma0::pi0selConfigurables pi0selOpts;
  o2::pwglf::sigma0::axisConfigurables axisOpts;

  o2::pwglf::sigma0::Sigma0BuilderModule sigma0Module;

  // For manual sliceBy  
  SliceCache cache;
  PresliceUnsorted<V0StandardDerivedDatas> perCollisionSTDDerived = o2::aod::v0data::straCollisionId;  

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> fverbose{"fverbose", true, "QA printout."};

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Initialize task
    sigma0Module.init(histos, evselOpts, lambdaselOpts, photonselOpts, sigma0selOpts, pi0selOpts, axisOpts);
  }

  // Dummy process function 
  void process(StraColls const&){} 

  void processRealData(StraColls const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&){
   auto start = std::chrono::high_resolution_clock::now();

   sigma0Module.process(collisions, fullV0s, histos, cache, ccdb, rateFetcher); 
   
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;

   if (fverbose) LOGF(info, "[Process function call, Sorted] N. Collisions: %i, N. V0s: %i, Processing time (s): %lf", collisions.size(), fullV0s.size(), elapsed.count());      
  }
  
  PROCESS_SWITCH(sigma0Unsorted, processRealData, "process run 3 real data", true);    
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
      adaptAnalysisTask<sigma0Sorted>(cfgc),
      adaptAnalysisTask<sigma0Unsorted>(cfgc)};
}
