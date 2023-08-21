// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseProducerMCTruthTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch

#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePhiSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Math/Vector4D.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "TMath.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{


using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;

// using FilteredFullV0s = soa::Filtered<aod::V0Datas>; /// predefined Join
// table for o2::aod::V0s = soa::Join<o2::aod::TransientV0s, o2::aod::StoredV0s>
// to be used when we add v0Filter
} // namespace o2::aod

/// \todo fix how to pass array to setSelection, getRow() passing a different
/// type!
// static constexpr float arrayV0Sel[3][3] = {{100.f, 100.f, 100.f}, {0.2f,
// 0.2f, 0.2f}, {100.f, 100.f, 100.f}}; unsigned int rows = sizeof(arrayV0Sel) /
// sizeof(arrayV0Sel[0]); unsigned int columns = sizeof(arrayV0Sel[0]) /
// sizeof(arrayV0Sel[0][0]);

template <typename T>
int getRowDaughters(int daughID, T const& vecID)
{
  int rowInPrimaryTrackTableDaugh = -1;
  for (size_t i = 0; i < vecID.size(); i++) {
    if (vecID.at(i) == daughID) {
      rowInPrimaryTrackTableDaugh = i;
      break;
    }
  }
  return rowInPrimaryTrackTableDaugh;
}

struct femtoUniverseProducerMCTruthTask {

int evCount = 0;
  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  /*
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
*/
  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};
  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};
  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> ConfEvtUseTPCmult{"ConfEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", true, "Activate filling of V0 into femtouniverse tables"};
  Configurable<bool> ConfIsActivatePhi{"ConfIsActivatePhi", true, "Activate filling of Phi into femtouniverse tables"};


  Configurable<std::vector<int>> ConfPDGCodes{"ConfPDGCodes", std::vector<int>{211,-211,2212,-2212}, "PDG of particles to be stored"};
  Configurable<bool> ConfAnalysisWithPID{"ConfAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles"};

  FemtoUniverseTrackSelection trackCuts;

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfUseGlobalTrack{"ConfUseGlobalTrack", true, "Use global track cuts (as defined in the official o2 documentation)"};
    Configurable<float> ConfPtLowGlobalTrack{"ConfPtLowGlobalTrack", 0.14, "Lower limit for Pt for the global track"};             // pT low
    Configurable<float> ConfPtHighGlobalTrack{"ConfPtHighGlobalTrack", 5.0, "Higher limit for Pt for the global track"};           // pT high
    Configurable<float> ConfEtaLowGlobalTrack{"ConfEtaLowGlobalTrack", -0.8, "Lower limit for Eta for the global track"};          // eta low
    Configurable<float> ConfEtaHighGlobalTrack{"ConfEtaHighGlobalTrack", 0.8, "Higher limit for Eta for the global track"};        // eta high
    Configurable<float> ConfDcaXYGlobalTrack{"ConfDcaXYGlobalTrack", 2.4, "Value for DCA_XY for the global track"};                // max dca to vertex XY
    Configurable<float> ConfDcaZGlobalTrack{"ConfDcaZGlobalTrack", 3.2, "Value for DCA_Z for the global track"};                   // max dca to vertex Z
    Configurable<int> ConfTpcClGlobalTrack{"ConfTpcClGlobalTrack", 88, "Number of tpc clasters for the global track"};             // min number of found TPC clusters
    Configurable<int> ConfTpcCrosRoGlobalTrack{"ConfTpcCrosRoGlobalTrack", 70, "Number of tpc crossed rows for the global track"}; // min number of crossed rows
    Configurable<float> ConfChi2TpcGlobalTrack{"ConfChi2TpcGlobalTrack", 4.0, "Chi2 / cluster for the TPC track segment for the global track"};
    Configurable<float> ConfChi2ItsGlobalTrack{"ConfChi2ItsGlobalTrack", 36.0, "Chi2 / cluster for the ITS track segment for the global track"};
  } ConfGlobalTracks;

   Configurable<std::vector<float>> ConfTrkCharge{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{
    FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{70.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.83f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kDCAzMax, "Track selection: ")}; /// \todo Reintegrate PID to the general selection container
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
  // Numbers from ~/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h //static constexpr ID Pion = 2; static constexpr ID Kaon = 3; static constexpr ID Proton = 4; static constexpr ID Deuteron = 5;


  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    if ((doprocessTrackMC) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }


    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);

    
    colCuts.init(&qaRegistry);
    
    //trackCuts.setSelection(ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    //trackCuts.setSelection(ConfTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    //trackCuts.setSelection(ConfTrkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    //trackCuts.setSelection(ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    /*
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

    
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoUniverseTrackSelection::kTPCfClsMin, femtoUniverseSelection::kLowerLimit);
    //--test
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoUniverseTrackSelection::kTPCcRowsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoUniverseTrackSelection::kTPCsClsMax, femtoUniverseSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femtoUniverseTrackSelection::kITSnClsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoUniverseTrackSelection::kITSnClsIbMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoUniverseTrackSelection::kDCAxyMax, femtoUniverseSelection::kAbsUpperLimit);
    // trackCuts.setSelection(ConfTrkCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoUniverseTrackSelection::kpTMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkPtmax, femtoUniverseTrackSelection::kpTMax, femtoUniverseSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkEta, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);

    //---end test
    trackCuts.setSelection(ConfTrkDCAzMax, femtoUniverseTrackSelection::kDCAzMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkPIDspecies);
    */
    trackCuts.setnSigmaPIDOffset(ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);



    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  // maybe in the future use track.passedPtRange() https://aliceo2group.github.io/analysis-framework/docs/basics-usage/HelperTasks.html?highlight=passedPtRange
  template <typename ParticleType>
  bool isGlobalTrack(ParticleType const& track)
  {
    if (track.pt() < ConfGlobalTracks.ConfPtLowGlobalTrack || track.pt() > ConfGlobalTracks.ConfPtHighGlobalTrack) { // pT cut
      return false;
    } else if (track.eta() < ConfGlobalTracks.ConfEtaLowGlobalTrack || track.eta() > ConfGlobalTracks.ConfEtaHighGlobalTrack) { // eta cut
      return false;
    } 
    return true;
  }




  template <typename CollisionType, typename TrackType>
  void fillCollisions(CollisionType const& col, TrackType const& tracks)
  {
   
    const auto vtxZ = col.posZ();
    const auto spher = 0; //colCuts.computeSphericity(col, tracks);
    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in
                                     /// FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfEvtUseTPCmult) {
      multNtr = col.multTPC();
    }
     // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions

    //CHECK WHAT CUTS SHOULD BE USED FOR MC TRUTH
    // if (!colCuts.isSelected(col)) {
    //   if (ConfIsTrigger) {
    //     outputCollision(vtxZ, mult, multNtr, spher, mMagField);
    //   }
    //   return;
    // }
    colCuts.fillQA(col);
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
  }

  template <typename TrackType>
  void fillParticles(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children

    for (auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if(!particle.isPhysicalPrimary()) continue;
      if(particle.eta()<-0.9 || particle.eta()>0.9) continue;
      if(particle.pt()<0.1 || particle.pt()>10) continue;

      uint32_t pdgCode = (uint32_t)particle.pdgCode();
      if(ConfAnalysisWithPID)
      {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfPDGCodes; // necessary due to some features of the Configurable
        for(uint32_t pdg: tmpPDGCodes){
          //LOGF(info,"%d %d",pdg,pdgCode);
          if((int)pdg==(int)pdgCode) pass=true;
        }
        if(!pass) continue;
      }
      //we cannot use isSelectedMinimal since it takes Ncls
      //if (!trackCuts.isSelectedMinimal(track)) {
      //  continue;
      //} 
 
      //trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
      //                 aod::femtouniverseparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      //auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(track);
      //instead of the bitmask, the PDG of the particle is stored as uint32_t

      

      // now the table is filled
       outputParts(outputCollision.lastIndex(), particle.pt(), particle.eta(),
                  particle.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                  0,
                  pdgCode,
                  0, childIDs, 0, 0);

    }
  }

  void
    processTrackMC(aod::FemtoFullCollisionMC const& col,
                   aod::BCsWithTimestamps const&,
                   aod::McCollisions const& mcCollisions,
                   aod::McParticles const& mcParticles)
  {
    if(evCount>30) return;
    evCount++;
    // magnetic field for run not needed for mc truth
    
    // fill the tables
    fillCollisions(col, mcParticles);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(femtoUniverseProducerMCTruthTask, processTrackMC, "Provide MC data for track analysis", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerMCTruthTask>(cfgc)};
  return workflow;
}

