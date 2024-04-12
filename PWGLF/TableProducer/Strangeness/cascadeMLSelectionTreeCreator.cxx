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
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Cascade ML Tree Creator task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//  This task loops over a set of V0 indices and
//  creates a TTree for ML training and testing.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//


#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct cascadeMLSelectionTreeCreator{
  Produces<aod::CascMLCandidates> cascadeMLCandidates;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // selection switches: please, use just one at a time! 
  Configurable<bool> getXiMinus{"getXiMinus", false, "If True, apply cuts to select XiMinus signal and bkg candidates"};
  Configurable<bool> getXiPlus{"getXiPlus", false, "If True, apply cuts to select XiPlus signal and bkg candidates"};
  Configurable<bool> getOmegaMinus{"getOmegaMinus", false, "If True, apply cuts to select OmegaMinus signal and bkg candidates"};
  Configurable<bool> getOmegaPlus{"getOmegaPlus", false, "If True, apply cuts to select OmegaPlus signal and bkg candidates"};

  Configurable<float> xiMassWindow{"xiMassWindow", 0.060, "Xi Mass Window around mass peak to consider (+/-, in GeV/c^2)"};
  Configurable<float> omegaMassWindow{"omegaMassWindow", 0.060, "Omega Mass Window around mass peak to consider (+/-, in GeV/c^2)"};

  // Master switch for selecting candidates
  struct : ConfigurableGroup {
    Configurable<float> Lambdav0cospa{"Lambdav0cospa", 0.90, "min V0 CosPA"};
    Configurable<float> Lambdadcav0dau{"Lambdadcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> Lambdadcanegtopv{"Lambdadcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> Lambdadcapostopv{"Lambdadcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> Lambdav0radius{"Lambdav0radius", 1.5, "minimum V0 radius (cm)"};
    Configurable<float> LambdaWindow{"LambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};
  } topoConfiguration;

  // Axis:
  ConfigurableAxis centralityAxis{"centralityAxis", {100, 0.0f, 100.0f}, ""};

  void init(InitContext const&)
  {
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {centralityAxis});
  }

    // Helper struct to pass v0 information
  struct {
    uint8_t massWindows;
    int charge;
    float centrality;

    // tracking properties
    int posITSCls;
    int negITSCls;
    int bachITSCls;
    uint32_t posITSClsSizes;
    uint32_t negITSClsSizes;
    uint32_t bachITSClsSizes;
    float posTPCRows;
    float negTPCRows;
    float bachTPCRows;
    
    // PID properties
    float posTPCSigmaPi;
    float negTPCSigmaPi;
    float posTPCSigmaPr;
    float negTPCSigmaPr;
    float bachTPCSigmaPi;
    float bachTPCSigmaKa;
    float TOFNSigmaXiLaPi;
    float TOFNSigmaXiLaPr;
    float TOFNSigmaXiPi;
    float TOFNSigmaOmLaPi;
    float TOFNSigmaOmLaPr;
    float TOFNSigmaOmKa;

    // Basic kine
    float mXi;
    float mOmega;
    float yXi;
    float yOmega;
    float mLambdaDaughter;
    float pt;
    float posEta;
    float negEta;
    float bachEta;

    // Topological
    float v0radius;
    float cascradius;
    float dcapostopv;
    float dcanegtopv;
    float dcabachtopv;
    float dcaV0Daughters;
    float dcaCascDaughters;
    float dcav0topv;
    float v0CosPA;
    float cascCosPA;

    // reserved for MC operation
    bool isXiMinus;
    bool isXiPlus;
    bool isOmegaMinus;
    bool isOmegaPlus; 
  } cascadeCandidate;
  
  // Process candidate and store properties in object
  template <typename TCollision, typename TCascObject>
  void processCandidate(TCollision const& coll, TCascObject const& cand)
  {    
    auto posTrackExtra = cand.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<dauTracks>();
    auto bachTrackExtra = cand.template bachTrackExtra_as<dauTracks>();

    cascadeCandidate.charge = cand.sign();
    cascadeCandidate.centrality = coll.centFT0C();

    // Track quality
    cascadeCandidate.posITSCls = posTrackExtra.itsNCls();
    cascadeCandidate.negITSCls = negTrackExtra.itsNCls();
    cascadeCandidate.bachITSCls = negTrackExtra.itsNCls();
    cascadeCandidate.posITSClsSizes = posTrackExtra.itsClusterSizes();
    cascadeCandidate.negITSClsSizes = negTrackExtra.itsClusterSizes();
    cascadeCandidate.bachITSClsSizes = negTrackExtra.itsClusterSizes();
    cascadeCandidate.posTPCRows = posTrackExtra.tpcCrossedRows();
    cascadeCandidate.negTPCRows = negTrackExtra.tpcCrossedRows();
    cascadeCandidate.bachTPCRows = bachTrackExtra.tpcCrossedRows();

    // TPC PID
    cascadeCandidate.posTPCSigmaPi = posTrackExtra.tpcNSigmaPi();
    cascadeCandidate.negTPCSigmaPi = negTrackExtra.tpcNSigmaPi();
    cascadeCandidate.posTPCSigmaPr = posTrackExtra.tpcNSigmaPr();
    cascadeCandidate.negTPCSigmaPr = negTrackExtra.tpcNSigmaPr();
    cascadeCandidate.bachTPCSigmaPi = bachTrackExtra.tpcNSigmaPi();
    cascadeCandidate.bachTPCSigmaKa = bachTrackExtra.tpcNSigmaKa();

    // TOF PID
    cascadeCandidate.TOFNSigmaXiLaPi = cand.tofNSigmaXiLaPi();
    cascadeCandidate.TOFNSigmaXiLaPr = cand.tofNSigmaXiLaPr();
    cascadeCandidate.TOFNSigmaXiPi = cand.tofNSigmaXiPi();
    cascadeCandidate.TOFNSigmaOmLaPi = cand.tofNSigmaOmLaPi();
    cascadeCandidate.TOFNSigmaOmLaPr = cand.tofNSigmaOmLaPr();
    cascadeCandidate.TOFNSigmaOmKa = cand.tofNSigmaOmKa();

    // Basic Kine
    cascadeCandidate.mXi = cand.mXi();
    cascadeCandidate.mOmega = cand.mOmega();
    cascadeCandidate.yXi = cand.yXi();
    cascadeCandidate.yOmega = cand.yOmega();
    cascadeCandidate.mLambdaDaughter = cand.mLambda();
    cascadeCandidate.pt = cand.pt();
    cascadeCandidate.posEta = cand.positiveeta();
    cascadeCandidate.negEta = cand.negativeeta();
    cascadeCandidate.bachEta = cand.bacheloreta();

    // Topological 
    cascadeCandidate.v0radius = cand.v0radius();
    cascadeCandidate.cascradius = cand.cascradius();
    cascadeCandidate.dcapostopv = cand.dcapostopv();
    cascadeCandidate.dcanegtopv = cand.dcanegtopv();
    cascadeCandidate.dcabachtopv = cand.dcabachtopv();
    cascadeCandidate.dcaV0Daughters = cand.dcaV0daughters();
    cascadeCandidate.dcaCascDaughters = cand.dcacascdaughters();
    cascadeCandidate.dcav0topv = cand.dcav0topv(coll.posX(), coll.posY(), coll.posZ());
    cascadeCandidate.v0CosPA = cand.v0cosPA(coll.posX(), coll.posY(), coll.posZ());
    cascadeCandidate.cascCosPA = cand.casccosPA(coll.posX(), coll.posY(), coll.posZ());

    // MC flags
    cascadeCandidate.isXiMinus = false;
    cascadeCandidate.isXiPlus = false;
    cascadeCandidate.isOmegaMinus = false;
    cascadeCandidate.isOmegaPlus = false;

    if constexpr ( requires {cand.pdgCode();} ) {
      cascadeCandidate.isXiMinus = (cand.pdgCode() == 3312);
      cascadeCandidate.isXiPlus = (cand.pdgCode() == -3312);
      cascadeCandidate.isOmegaMinus = (cand.pdgCode() == 3334);
      cascadeCandidate.isOmegaPlus = (cand.pdgCode() == -3334);
    }

    // Mass window selections
    cascadeCandidate.massWindows = 0;
    if(std::abs(cascadeCandidate.mXi-1.322) < xiMassWindow) {
      cascadeCandidate.massWindows = cascadeCandidate.massWindows | 1 << 0 ; 
    }
    if(std::abs(cascadeCandidate.mOmega-1.6725) < omegaMassWindow) {
      cascadeCandidate.massWindows = cascadeCandidate.massWindows | 1 << 1 ; 
    }
    if(cascadeCandidate.massWindows == 0){
      return; // skip
    }

    // populate
    cascadeMLCandidates(
      cascadeCandidate.massWindows,
      cascadeCandidate.charge,
      cascadeCandidate.centrality,

      // Track quality
      cascadeCandidate.posITSCls,
      cascadeCandidate.negITSCls,
      cascadeCandidate.bachITSCls,
      cascadeCandidate.posITSClsSizes,
      cascadeCandidate.negITSClsSizes,
      cascadeCandidate.bachITSClsSizes,
      cascadeCandidate.posTPCRows,
      cascadeCandidate.negTPCRows,
      cascadeCandidate.bachTPCRows,

      // TPC PID
      cascadeCandidate.posTPCSigmaPi,
      cascadeCandidate.negTPCSigmaPi,
      cascadeCandidate.posTPCSigmaPr,
      cascadeCandidate.negTPCSigmaPr,
      cascadeCandidate.bachTPCSigmaPi,
      cascadeCandidate.bachTPCSigmaKa,

      // TOF PID
      cascadeCandidate.TOFNSigmaXiLaPi,
      cascadeCandidate.TOFNSigmaXiLaPr,
      cascadeCandidate.TOFNSigmaXiPi,
      cascadeCandidate.TOFNSigmaOmLaPi,
      cascadeCandidate.TOFNSigmaOmLaPr,
      cascadeCandidate.TOFNSigmaOmKa,

      // Basic Kine
      cascadeCandidate.mXi,
      cascadeCandidate.mOmega,
      cascadeCandidate.yXi,
      cascadeCandidate.yOmega,
      cascadeCandidate.mLambdaDaughter,
      cascadeCandidate.pt,
      cascadeCandidate.posEta,
      cascadeCandidate.negEta,
      cascadeCandidate.bachEta,

      // Topological 
      cascadeCandidate.v0radius,
      cascadeCandidate.cascradius,
      cascadeCandidate.dcapostopv,
      cascadeCandidate.dcanegtopv,
      cascadeCandidate.dcabachtopv,
      cascadeCandidate.dcaV0Daughters,
      cascadeCandidate.dcaCascDaughters,
      cascadeCandidate.dcav0topv,
      cascadeCandidate.v0CosPA,
      cascadeCandidate.cascCosPA,

      // MC identity
      cascadeCandidate.isXiMinus,
      cascadeCandidate.isXiPlus,
      cascadeCandidate.isOmegaMinus,
      cascadeCandidate.isOmegaPlus
    );
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascTOFNSigmas> const& cascades, dauTracks const&)
  {
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    for (auto& casc: cascades){ // looping over lambdas 
      processCandidate(coll, casc);
    }
  }
  void processSimData(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascMCDatas, aod::CascTOFNSigmas> const& cascades, dauTracks const&)
  {
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    for (auto& casc: cascades){ // looping over lambdas 
      processCandidate(coll, casc);
    }
  }

  PROCESS_SWITCH(cascadeMLSelectionTreeCreator, processRealData, "Produce Run 3 cascade tables (real data)", false);  
  PROCESS_SWITCH(cascadeMLSelectionTreeCreator, processSimData, "Produce Run 3 cascade tables (simulation)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cascadeMLSelectionTreeCreator>(cfgc)};
  
}