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

/// \file femtoDreamProducerReducedTask.cxx
/// \brief Tasks that produces the track tables used for the pairing (tracks only)
/// \author Luca Barioglio, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <vector>
#include "TMath.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Math/Vector4D.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
using FemtoFullCollision_noCent_MC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;

using FemtoFullTracks = soa::Join<aod::FullTracks, aod::TracksDCA,
                                  aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa,
                                  aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                                  aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa,
                                  aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;
} // namespace o2::aod

struct femtoDreamProducerReducedTask {

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or Run2"}; // Choose if running on converted data or Run3
  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  // Event cuts
  FemtoDreamCollisionSelection colCuts;
  Configurable<bool> ConfEvtUseTPCmult{"ConfEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfEvtAddOfflineCheck{"ConfEvtAddOfflineCheck", false, "Evt sel: additional checks for offline selection (not part of sel8 yet)"};

  Configurable<bool> ConfTrkRejectNotPropagated{"ConfTrkRejectNotPropagated", false, "True: reject not propagated tracks"};

  Configurable<int> ConfTrkPDGCode{"ConfTrkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  // Track cuts
  FemtoDreamTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")}; /// here we need an open cut to do the DCA fits later on!
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "Conf"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  // off set the center of the nsigma distribution to deal with bad TPC/TOF calibration
  Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry Registry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {

    int CutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    Registry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});

    colCuts.setCuts(ConfEvtZvtx.value, ConfEvtTriggerCheck.value, ConfEvtTriggerSel.value, ConfEvtOfflineCheck.value, ConfEvtAddOfflineCheck.value, ConfIsRun3.value);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack,
                   aod::femtodreamparticle::TrackType::kNoChild,
                   aod::femtodreamparticle::cutContainerType>(&qaRegistry, &Registry);
    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // get magnetic field for run
    if (mRunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();
    float output = -999;

    if (ConfIsRun3 && !ConfIsForceGRP) {
      static o2::parameters::GRPMagField* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
      // taken from GRP onject definition of getNominalL3Field; update later to something smarter (mNominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      output = 0.1 * (NominalL3Field);

    } else {

      static o2::parameters::GRPObject* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
      output = 0.1 * (grpo->getNominalL3Field());
    }
    mMagField = output;
    mRunNumber = bc.runNumber();
  }

  template <typename CollisionType, typename CandidateType>
  void fillMCParticle(CollisionType const& col, CandidateType const& candidate, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (candidate.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = candidate.mcParticle();
      auto pdgCode = particleMC.pdgCode();

      int particleOrigin = 99;
      auto motherparticleMC = particleMC.template mothers_as<aod::McParticles>().front();

      if (abs(pdgCode) == abs(ConfTrkPDGCode.value)) {
        if ((col.has_mcCollision() && (particleMC.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
        } else if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (motherparticleMC.isPhysicalPrimary() && particleMC.getProcess() == 4) {
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
        } else if (particleMC.getGenStatusCode() == -1) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
      } else {

        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
    } else {
      outputPartsMCLabels(-1);
    }
  }

  template <bool isMC, bool useCentrality, typename CollisionType, typename TrackType>
  void fillCollisionsAndTracks(CollisionType const& col, TrackType const& tracks) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      if constexpr (useCentrality) {
        mult = col.centFT0M();
      } else {
        mult = 0.;
      }
      multNtr = col.multNTracksPV();
    } else {
      mult = 1; // multiplicity percentile is known in Run 2
      multNtr = col.multTracklets();
    }
    if (ConfEvtUseTPCmult) {
      multNtr = col.multTPC();
    }
    colCuts.fillQA(col, mult);

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    /// That includes checking if there are any usable tracks in a collision
    if (!colCuts.isSelectedCollision(col)) {
      return;
    }
    if (colCuts.isEmptyCollision(col, tracks, trackCuts)) {
      return;
    }

    // now the table is filled
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);

    // these IDs are necessary to keep track of the children
    // since this producer only produces the tables for tracks, there are no children
    std::vector<int> childIDs = {0, 0};
    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild>(track);
      // an array of two bit-wise containers of the systematic variations is obtained
      // one container for the track quality cuts and one for the PID cuts
      auto cutContainer = trackCuts.getCutContainer<true, aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtodreamparticle::ParticleType::kTrack,
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);
      if constexpr (isMC) {
        fillMCParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }

      if (ConfIsDebug) {
        outputDebugParts(track.sign(),
                         (uint8_t)track.tpcNClsFound(),
                         track.tpcNClsFindable(),
                         (uint8_t)track.tpcNClsCrossedRows(),
                         track.tpcNClsShared(),
                         track.tpcInnerParam(),
                         track.itsNCls(),
                         track.itsNClsInnerBarrel(),
                         track.dcaXY(),
                         track.dcaZ(),
                         track.tpcSignal(),
                         track.tpcNSigmaEl(),
                         track.tpcNSigmaPi(),
                         track.tpcNSigmaKa(),
                         track.tpcNSigmaPr(),
                         track.tpcNSigmaDe(),
                         track.tpcNSigmaTr(),
                         track.tpcNSigmaHe(),
                         track.tofNSigmaEl(),
                         track.tofNSigmaPi(),
                         track.tofNSigmaKa(),
                         track.tofNSigmaPr(),
                         track.tofNSigmaDe(),
                         track.tofNSigmaTr(),
                         track.tofNSigmaHe(),
                         -1,
                         track.itsNSigmaEl(),
                         track.itsNSigmaPi(),
                         track.itsNSigmaKa(),
                         track.itsNSigmaPr(),
                         track.itsNSigmaDe(),
                         track.itsNSigmaTr(),
                         track.itsNSigmaHe(),
                         -999., -999., -999., -999., -999., -999.);
      }
    }
  }

  void
    processData(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);
    // fill the tables
    fillCollisionsAndTracks<false, true>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(femtoDreamProducerReducedTask, processData, "Provide experimental data", true);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::McCollisions const&, aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    auto tracksWithItsPid = soa::Attach<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);
    // fill the tables
    fillCollisionsAndTracks<true, true>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(femtoDreamProducerReducedTask, processMC, "Provide MC data", false);

  void processMC_noCentrality(aod::FemtoFullCollision_noCent_MC const& col,
                              aod::BCsWithTimestamps const&,
                              soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                              aod::McCollisions const&, aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    auto tracksWithItsPid = soa::Attach<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);
    // fill the tables
    fillCollisionsAndTracks<true, false>(col, tracksWithItsPid);
  }
  PROCESS_SWITCH(femtoDreamProducerReducedTask, processMC_noCentrality, "Provide MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamProducerReducedTask>(cfgc)};
  return workflow;
}
