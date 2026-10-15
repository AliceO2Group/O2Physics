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
/// \file dielptonEffi.cxx
/// \brief calculating dilepton efficiency for dalitz analyses
/// \author Laura Gans-Bartl, laura.gans-bartl@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include <PWGEM/Dilepton/Utils/MCUtilities.h>

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;

using MyEMPrimaryElectrons = soa::Join<aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronMCLabels, aod::EMPrimaryElectronDaEMEventIds_001>;
using MyEMPrimaryElectron = MyEMPrimaryElectrons::iterator;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::pidTOFbeta, aod::TOFSignal, aod::TOFEvTime>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels, aod::mcTPCTuneOnData>;
using MyTrackMC = MyTracksMC::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyCollisionsMC = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels, o2::aod::EmMagFields>;
using MyCollisionMC = MyCollisionsMC::iterator;

struct DileptonEffiForDalitzSkimmer {
  SliceCache cache;
  PresliceOptional<MyEMPrimaryElectrons> perCollision = aod::emprimaryelectronda::pmeventId;
  PresliceOptional<MyTracksMC> perTracksCollision = aod::track::collisionId;
  PresliceOptional<aod::McParticles> perDileptonMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0KFEMEventIds> perCollisionV0 = aod::v0photonkf::pmeventId;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  bool isAmbiguousTrack(const auto& track, const aod::AmbiguousTracks& ambiguousTracks)
  {
    for (const auto& ambiguousTrack : ambiguousTracks) {
      if (ambiguousTrack.trackId() == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  void init(o2::framework::InitContext&)
  {
    fRegistry.add("hPtSglEpEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEmEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEpEMPrimTruePtRec", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEmEMPrimTruePtRec", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSglPipEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglPimEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSglEpMCPart", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEmMCPart", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSgltpcNSigmaEl", "pT;p_{T} (GeV/c); tpcNSigmaEl", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);
    fRegistry.add("hPtSgltpcNSigmaPi", "pT;p_{T} (GeV/c); tpcNSigmaPi", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);
    fRegistry.add("hPtSgltofNSigmaEl", "pT;p_{T} (GeV/c); tofNSigmaEl", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);

    fRegistry.add("hPtSglTracktpcNSigmaEl", "pT;p_{T} (GeV/c); tpcNSigmaEl", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);
    fRegistry.add("hPtSglTracktpcNSigmaPi", "pT;p_{T} (GeV/c); tpcNSigmaPi", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);
    fRegistry.add("hPtSglTracktofNSigmaEl", "pT;p_{T} (GeV/c); tofNSigmaEl", kTH2F, {{1000, 0.0f, 10}, {200, -5.f, 5}}, false);

    fRegistry.add("hPtMcEtaGeneral", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtMcEtaDalitz", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtEPrimMulti", ";p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hEPrimNumOccurence", ";p_{T} (GeV/c)", kTH1F, {{11, -0.5f, 10.5}}, false);

    fRegistry.add("hPtTrackMulti", ";p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hTrackNumOccurence", ";p_{T} (GeV/c)", kTH1F, {{11, -0.5f, 10.5}}, false);

    fRegistry.add("hPtSglEpTrackEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEmTrackEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSglEmFromTrueEtaGen", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtGammaFromTrueEtaGen", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSglEmFromTrueEtaRecoRecoPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglEmFromTrueEtaRecoGenPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtGammaFromTrueEtaReco", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hPtSglPipTrackEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("hPtSglPimTrackEMPrimTrue", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);

    fRegistry.add("hMcPartsPerColl", " ;N_{tracks}/collision; N", kTH1F, {{400, -0.5, 400.5}}, false);

    fRegistry.add("hCheckAsso", "pT;p_{T} (GeV/c)", kTH1F, {{2, -0.5f, 1.5}}, false);
  }
  void processReco(aod::PMEvents const& collisions, MyEMPrimaryElectrons const& emprimaryelectrons, aod::EMMCParticles_001 const& emmcparticle) // aod::EMEvents_005, aod::EMMCParticles_001
  {
    LOG(info) << "started process";
    if (emprimaryelectrons.size() <= 0 || emmcparticle.size() <= 0) {
      return;
    }
    LOG(info) << "size emprimaryelectrons: " << emprimaryelectrons.size();

    const int pdgEp = -11;
    const int pdgEm = 11;
    const int pdgPip = -211;
    const int pdgPim = 211;
    auto iterMC = emmcparticle.begin();
    std::unordered_map<uint64_t, uint32_t> mapMultCount;

    //**Single electrons and positrons**/
    for (auto const& collision : collisions) {
      auto slicedParticles = emprimaryelectrons.sliceBy(perCollision, collision.globalIndex());

      for (const auto& ePrimary : slicedParticles) {
        if (!std::isfinite(ePrimary.pt())) { // || std::abs(ePrimary.eta()) > 0.8
          continue;
        }

        fRegistry.fill(HIST("hPtSgltpcNSigmaEl"), ePrimary.pt(), ePrimary.tpcNSigmaEl());
        fRegistry.fill(HIST("hPtSgltpcNSigmaPi"), ePrimary.pt(), ePrimary.tpcNSigmaPi());
        fRegistry.fill(HIST("hPtSgltofNSigmaEl"), ePrimary.pt(), ePrimary.tofNSigmaEl());

        const auto mcId = ePrimary.emmcparticleId();

        if (mcId < 0 || mcId >= emmcparticle.size()) {
          continue;
        }
        iterMC.setCursor(mcId);
        mapMultCount[static_cast<uint64_t>(iterMC.globalIndex())]++;

        if (!iterMC.isPhysicalPrimary() || fabs(iterMC.y()) > 0.8) {
          continue;
        }

        if (iterMC.pdgCode() == pdgEp) {
          fRegistry.fill(HIST("hPtSglEpEMPrimTrue"), iterMC.pt());
          fRegistry.fill(HIST("hPtSglEpEMPrimTruePtRec"), ePrimary.pt());
        }
        if (iterMC.pdgCode() == pdgEm) {
          fRegistry.fill(HIST("hPtSglEmEMPrimTrue"), iterMC.pt());
          fRegistry.fill(HIST("hPtSglEmEMPrimTruePtRec"), ePrimary.pt());
        }
        if (iterMC.pdgCode() == pdgPip) {
          fRegistry.fill(HIST("hPtSglPipEMPrimTrue"), iterMC.pt());
        }
        if (iterMC.pdgCode() == pdgPim) {
          fRegistry.fill(HIST("hPtSglPimEMPrimTrue"), iterMC.pt());
        }
      }
    }
    uint32_t maxNumOccurrences = 0;
    uint32_t mcIndexWithMaxTracks = 0;
    uint32_t numOkMc = 0;
    uint32_t numNotOkMc = 0;

    for (const auto& [mcGlobalIndex, numOccurrences] : mapMultCount) {
      if (numOccurrences > 1) {
        iterMC.setCursor(mcGlobalIndex);
        fRegistry.fill(HIST("hPtEPrimMulti"), iterMC.pt());
        fRegistry.fill(HIST("hEPrimNumOccurence"), numOccurrences);
        numNotOkMc++;
      } else {
        numOkMc++;
      }
      if (numOccurrences > maxNumOccurrences) {
        maxNumOccurrences = numOccurrences;
        mcIndexWithMaxTracks = mcGlobalIndex;
      }
    }
    LOG(info) << "numOkMc: " << numOkMc;
    LOG(info) << "numNotOkMc: " << numNotOkMc;
  }
  PROCESS_SWITCH(DileptonEffiForDalitzSkimmer, processReco, "process for reco table emprimary", false);

  void processTracks(aod::Collisions const& collisions, aod::McParticles const& mcParticles, MyTracksMC const& tracks)
  {
    if (tracks.size() <= 0 || mcParticles.size() <= 0) {
      return;
    }
    LOG(info) << "size tracks: " << tracks.size();

    ////**Analysing tracks **/
    const int pdgEp = -11;
    const int pdgEm = 11;
    const int pdgPip = -211;
    const int pdgPim = 211;
    auto iterMCTrack = mcParticles.begin();

    std::unordered_map<uint64_t, uint32_t> mapMultCount;
    ////******Tracks before process function******////
    for (auto const& collision : collisions) {
      auto slicedTracks = tracks.sliceBy(perTracksCollision, collision.globalIndex());
      for (const auto& track : slicedTracks) {
        if (!track.has_mcParticle() || !std::isfinite(track.pt())) { //|| std::abs(track.eta()) > 0.8
          continue;
        }

        fRegistry.fill(HIST("hPtSglTracktpcNSigmaEl"), track.pt(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("hPtSglTracktpcNSigmaPi"), track.pt(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("hPtSglTracktofNSigmaEl"), track.pt(), track.tofNSigmaEl());

        const auto mcId = track.mcParticleId();

        if (mcId < 0 || mcId >= mcParticles.size()) {
          continue;
        }

        iterMCTrack.setCursor(mcId);
        mapMultCount[static_cast<uint64_t>(iterMCTrack.globalIndex())]++;

        if (!iterMCTrack.isPhysicalPrimary() || fabs(iterMCTrack.y()) > 0.8) {
          continue;
        }

        if (iterMCTrack.pdgCode() == pdgEp) {
          fRegistry.fill(HIST("hPtSglEpTrackEMPrimTrue"), iterMCTrack.pt());
        }

        if (iterMCTrack.pdgCode() == pdgEm) {
          fRegistry.fill(HIST("hPtSglEmTrackEMPrimTrue"), iterMCTrack.pt());
        }
        if (iterMCTrack.pdgCode() == pdgPip) {
          fRegistry.fill(HIST("hPtSglPipTrackEMPrimTrue"), iterMCTrack.pt());
        }
        if (iterMCTrack.pdgCode() == pdgPim) {
          fRegistry.fill(HIST("hPtSglPimTrackEMPrimTrue"), iterMCTrack.pt());
        }
      }
    }
    uint32_t maxNumOccurrences = 0;
    uint64_t mcIndexWithMaxTracks = 0;
    uint32_t numOkMc = 0;
    uint64_t numNotOkMc = 0;

    for (const auto& [mcGlobalIndex, numOccurrences] : mapMultCount) {
      if (numOccurrences > 1) {
        iterMCTrack.setCursor(mcGlobalIndex);
        // LOG(info) << "Track | Global Index: " << iterMCTrack.globalIndex();
        fRegistry.fill(HIST("hPtTrackMulti"), iterMCTrack.pt());
        fRegistry.fill(HIST("hTrackNumOccurence"), numOccurrences);
        numNotOkMc++;
      } else {
        numOkMc++;
      }
      if (numOccurrences > maxNumOccurrences) {
        maxNumOccurrences = numOccurrences;
        mcIndexWithMaxTracks = mcGlobalIndex;
      }
    }
  }
  PROCESS_SWITCH(DileptonEffiForDalitzSkimmer, processTracks, "process for reco MyTracks", true);

  using FilteredMyCollisionsMC = soa::Filtered<MyCollisionsMC>;
  void processGens(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, aod::PMEvents const& pmcollisions, MyEMPrimaryElectrons const& emprimaryelectrons, aod::EMMCParticles_001 const& emmcparticle)
  {
    LOG(info) << "started process MC";
    if (mcParticles.size() <= 0 || emprimaryelectrons.size() <= 0 || emmcparticle.size() <= 0) {
      return;
    }
    LOG(info) << " size mcParticles: " << mcParticles.size();
    auto iterChild = mcParticles.begin();
    auto iterMom = mcParticles.begin();

    const int pdgEp = -11;
    const int pdgEm = 11;
    const int pdgEta = 221;
    const int pdgGamma = 22;

    ////******generated Dileptons******////
    std::vector<int64_t> vecEpMCParticles;
    vecEpMCParticles.reserve(1e5);
    std::vector<int64_t> vecEmMCParticles;
    vecEpMCParticles.reserve(1e5);

    std::vector<float> vecEmMCPartPt;
    vecEmMCPartPt.reserve(1e5);

    //**Single electrons and positrons**//
    for (auto const& mcCollision : mcCollisions) {
      auto slicedMcParticles = mcParticles.sliceBy(perDileptonMcCollision, mcCollision.globalIndex());
      fRegistry.fill(HIST("hMcPartsPerColl"), slicedMcParticles.size());
      for (const auto& mcParticle : slicedMcParticles) {
        if (!mcParticle.isPhysicalPrimary() || fabs(mcParticle.y()) > 0.8) {
          continue;
        }
        vecEmMCPartPt.push_back(mcParticle.pt());

        // const float dx = mcParticle.vx();  //- mcCollision.posX()
        // const float dy = mcParticle.vy();
        // const float dz = mcParticle.vz();

        // const float prodDist = sqrt(dx * dx + dy * dy);     //+ dz * dz

        if (mcParticle.pdgCode() == pdgEp) {
          if (!std::isfinite(mcParticle.pt())) {
            continue;
          }
          fRegistry.fill(HIST("hPtSglEpMCPart"), mcParticle.pt());
          vecEpMCParticles.push_back(mcParticle.globalIndex());
        }
        if (mcParticle.pdgCode() == pdgEm) {
          if (!std::isfinite(mcParticle.pt())) {
            continue;
          }
          fRegistry.fill(HIST("hPtSglEmMCPart"), mcParticle.pt());
          vecEmMCParticles.push_back(mcParticle.globalIndex());
        }

        ////******Investigating Eta Daltz Decay******////
        if (!mcParticle.has_mothers()) {
          continue;
        }
        if (mcParticle.mothersIds()[0] < 0 || mcParticle.mothersIds()[0] >= mcParticles.size()) {
          continue;
        }
        iterMom.setCursor(mcParticle.mothersIds()[0]);
        if (iterMom.pdgCode() != pdgEta) {
          continue;
        }
        fRegistry.fill(HIST("hPtMcEtaGeneral"), iterMom.pt());

        if (!iterMom.has_daughters()) {
          continue;
        }
        bool childEp = false, childEm = false, childG = false;
        double ptEmFromDalitzEta = -1;
        double ptGammaFromDalitzEta = -1;
        for (int childId = iterMom.daughtersIds()[0]; childId <= iterMom.daughtersIds()[1]; ++childId) {
          iterChild.setCursor(childId);
          if (iterChild.pdgCode() == pdgEp) {
            childEp = true;
          }
          if (iterChild.pdgCode() == pdgEm) {
            childEm = true;
            ptEmFromDalitzEta = iterChild.pt();
          }
          if (iterChild.pdgCode() == pdgGamma) {
            childG = true;
            ptGammaFromDalitzEta = iterChild.pt();
          }
        }
        if (childEp && childEm && childG) {
          fRegistry.fill(HIST("hPtMcEtaDalitz"), iterMom.pt());
          fRegistry.fill(HIST("hPtSglEmFromTrueEtaGen"), ptEmFromDalitzEta);
          fRegistry.fill(HIST("hPtGammaFromTrueEtaGen"), ptGammaFromDalitzEta);
        }
      }
      if (vecEpMCParticles.empty() || vecEmMCParticles.empty()) {
        continue;
      }
    }

    // Reco loop Elektronen
    auto iterMC = emmcparticle.begin();
    auto iterMomE = emmcparticle.begin();
    for (auto const& pmcollision : pmcollisions) {
      auto slicedParticles = emprimaryelectrons.sliceBy(perCollision, pmcollision.globalIndex());

      for (const auto& ePrimary : slicedParticles) {
        if (!std::isfinite(ePrimary.pt()) || fabs(iterMC.y()) > 0.8) { // || std::abs(ePrimary.eta()) > 0.8
          continue;
        }
        const auto mcId = ePrimary.emmcparticleId();
        iterMC.setCursor(mcId);

        if (!iterMC.isPhysicalPrimary() || fabs(iterMC.y()) > 0.8) {
          continue;
        }

        if (iterMC.pdgCode() != pdgEm) {
          continue;
        }

        if (!iterMC.has_mothers()) {
          continue;
        }
        iterMomE.setCursor(iterMC.mothersIds()[0]);
        if (iterMomE.pdgCode() != pdgEta) {
          continue;
        }
        fRegistry.fill(HIST("hPtSglEmFromTrueEtaRecoRecoPt"), ePrimary.pt());
        fRegistry.fill(HIST("hPtSglEmFromTrueEtaRecoGenPt"), iterMC.pt()); // to do!!!!!!!!!!!!!!!
      }
    }
  }
  PROCESS_SWITCH(DileptonEffiForDalitzSkimmer, processGens, "process for generated (MCParticles)", false);

  void processCheckAssociation(MyEMPrimaryElectrons const& emprimaryelectrons, MyTracksMC const& tracks)
  {
    auto iterTrack = tracks.begin();
    LOG(info) << "start process Association";
    for (const auto& ePrimary : emprimaryelectrons) {
      int idEm = ePrimary.emmcparticleId();
      iterTrack.setCursor(ePrimary.trackId());
      int idTrack = iterTrack.mcParticleId();

      if (idEm == idTrack) {
        fRegistry.fill(HIST("hCheckAsso"), 1);
      } else {
        fRegistry.fill(HIST("hCheckAsso"), 0);
      }
    }
  }
  PROCESS_SWITCH(DileptonEffiForDalitzSkimmer, processCheckAssociation, "checks if McParticle association is valid", false);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<DileptonEffiForDalitzSkimmer>(context, TaskName{"dilepton-effi-for-dalitz-skimmer"})};
}
