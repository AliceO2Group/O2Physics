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

/// \file femtoDreamProducer.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Ravindra Singh, GSI, ravindra.singh@cern.ch

#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct femtoDreamProducer {

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDHfCand> rowCandCharmHad;
  Produces<aod::FDHfCandMC> rowCandMCCharmHad;
  Produces<aod::FDHfCandMCGen> rowCandCharmHadGen;
  Produces<aod::FDParticlesIndex> outputPartsIndex;
  Produces<aod::FDMCCollisions> outputMCCollision;
  Produces<aod::FDMCCollLabels> outputCollsMCLabels;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
  Produces<aod::FDExtMCLabels> outputPartsExtMCLabels;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Configurable<bool> isDebug{"isDebug", true, "Enable Debug tables"};
  Configurable<bool> isRun3{"isRun3", true, "Running on Run3 or pilot"};
  Configurable<bool> isForceGRP{"isForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  /// Event selection
  Configurable<bool> doPvRefit{"doPvRefit", false, "do PV refit excluding the considered track"};
  Configurable<float> evtZvtx{"evtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> evtTriggerCheck{"evtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> evtTriggerSel{"evtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> evtOfflineCheck{"evtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> evtAddOfflineCheck{"evtAddOfflineCheck", false, "Evt sel: additional checks for offline selection (not part of sel8 yet)"};

  /// Lc table
  Configurable<bool> useCent{"useCent", false, "Enable centrality for lc"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};

  Configurable<bool> trkRejectNotPropagated{"trkRejectNotPropagated", false, "True: reject not propagated tracks"};
  Configurable<int> trkPDGCode{"trkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  Configurable<std::vector<float>> trkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "trk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> trkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "trk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> trkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "trk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> trkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "trk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "trk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "trk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "trk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "trk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "trk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "trk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "trk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "trk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> trkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "trk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> trkPIDnSigmaOffsetTPC{"trkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> trkPIDnSigmaOffsetTOF{"trkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton}, "Trk sel: Particles species for PID"};

  using CandidateLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;
  using CandidateLcMC = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>;

  using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
  using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
  using FemtoFullMCgenCollisions = soa::Join<aod::McCollisions, o2::aod::MultsExtraMC>;
  using FemtoFullMCgenCollision = FemtoFullMCgenCollisions::iterator;
  using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FemtoHFTrack = FemtoHFTracks::iterator;

  using GeneratedMC = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;

  FemtoDreamCollisionSelection colCuts;
  FemtoDreamTrackSelection trackCuts;

  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry TrackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};

  HfHelper hfHelper;
  int runNumber;
  float magField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  o2::base::MatLayerCylSet* lut;
  // if (doPvRefit){ lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));} //! may be it useful, will check later

  void init(InitContext&)
  {
    std::array<bool, 5> processes = {doprocessDataCharmHad, doprocessMCCharmHad, doprocessDataCharmHadWithML, doprocessMCCharmHadWithML, doprocessMCCharmHadGen};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int CutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    TrackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});

    colCuts.setCuts(evtZvtx.value, evtTriggerCheck.value, evtTriggerSel.value, evtOfflineCheck.value, evtAddOfflineCheck.value, isRun3.value);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(trkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(trkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(trkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(trkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(trkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(trkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(trkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(trkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(trkPIDspecies);
    trackCuts.setnSigmaPIDOffset(trkPIDnSigmaOffsetTPC, trkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &TrackRegistry);

    runNumber = 0;
    magField = 0.0;
    /// Initializing CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // auto bc = col.bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, !isRun3 ? ccdbPathGrp : ccdbPathGrpMag, lut, !isRun3);

    //    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    //    // get magnetic field for run
    //    if (runNumber == bc.runNumber())
    //      return;
    //    auto timestamp = bc.timestamp();
    //    float output = -999;
    //
    //    if (isRun3 && !isForceGRP) {
    //      static o2::parameters::GRPMagField* grpo = nullptr;
    //      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    //      if (grpo == nullptr) {
    //        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
    //        return;
    //      }
    //      LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
    //      // taken from GRP onject definition of getNominalL3Field; update later to something smarter (mNominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
    //      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
    //      output = 0.1 * (NominalL3Field);
    //
    //    } else {
    //
    //      static o2::parameters::GRPObject* grpo = nullptr;
    //      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
    //      if (grpo == nullptr) {
    //        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
    //        return;
    //      }
    //      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    //      output = 0.1 * (grpo->getNominalL3Field());
    //    }
    //    magField = output;
    //    runNumber = bc.runNumber();
  }

  template <typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    outputDebugParts(particle.sign(),
                     (uint8_t)particle.tpcNClsFound(),
                     particle.tpcNClsFindable(),
                     (uint8_t)particle.tpcNClsCrossedRows(),
                     particle.tpcNClsShared(),
                     particle.tpcInnerParam(),
                     particle.itsNCls(),
                     particle.itsNClsInnerBarrel(),
                     particle.dcaXY(),
                     particle.dcaZ(),
                     particle.tpcSignal(),
                     particle.tpcNSigmaPi(),
                     particle.tpcNSigmaKa(),
                     particle.tpcNSigmaPr(),
                     particle.tofNSigmaPi(),
                     particle.tofNSigmaKa(),
                     particle.tofNSigmaPr(),
                     -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.);
  }

  template <typename CollisionType, typename ParticleType>
  void fillMCParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();
      int particleOrigin = 99;
      int pdgCodeMother = -1;
      // get list of mothers
      // could be empty (for example in case of injected light nuclei)
      auto motherparticlesMC = particleMC.template mothers_as<aod::McParticles>();
      // check pdg code
      if (std::abs(pdgCode) == abs(trkPDGCode.value)) {
        if ((col.has_mcCollision() && (particleMC.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
        } else if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (particleMC.getGenStatusCode() == -1) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
        } else if (!motherparticlesMC.empty()) {
          // get direct mother of the particle
          auto motherparticleMC = motherparticlesMC.front();
          pdgCodeMother = motherparticleMC.pdgCode();
          if (motherparticleMC.isPhysicalPrimary() && particleMC.getProcess() == 4) {
            particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
          }
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
      } else {
        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
      if (isDebug) {
        outputPartsExtMCLabels(outputPartsMC.lastIndex());
        outputDebugPartsMC(pdgCodeMother);
      }
    } else {
      outputPartsMCLabels(-1);
      if (isDebug) {
        outputPartsExtMCLabels(-1);
      }
    }
  }

  template <typename CollisionType>
  void fillMCCollision(CollisionType const& col)
  {
    if (col.has_mcCollision()) {
      auto genMCcol = col.template mcCollision_as<FemtoFullMCgenCollisions>();
      outputMCCollision(genMCcol.multMCNParticlesEta08());
      outputCollsMCLabels(outputMCCollision.lastIndex());
    } else {
      outputCollsMCLabels(-1);
    }
  }

  template <bool isMC, typename TrackType>
  bool fillTracksForCharmHadron(TrackType const& tracks, FemtoHFTrack const& prong0, FemtoHFTrack const& prong1, FemtoHFTrack const& prong2, int candSize)
  {

    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    // std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    bool fIsTrackFilled = false;

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      if ((candSize == 1) && (track.globalIndex() == prong0.globalIndex() || track.globalIndex() == prong1.globalIndex() || track.globalIndex() == prong2.globalIndex()))
        continue;

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, true>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

      // track global index
      outputPartsIndex(track.globalIndex());
      // now the table is filled

      outputParts(outputCollision.lastIndex() + 1,
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtodreamparticle::ParticleType::kTrack,
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);
      fIsTrackFilled = true;
      // tmpIDtrack.push_back(track.globalIndex());
      if (isDebug.value) {
        fillDebugParticle(track);
      }

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtodreamparticle::ParticleType::kTrack);
        //    if constexpr (isMC) {
        //      fillMCParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
        //    }
      }
    }
    return fIsTrackFilled;
  }

  template <bool isMC, bool useCharmMl, typename TrackType, typename CollisionType, typename CandType>
  void fillCharmHadronTable(CollisionType const& col, TrackType const& tracks, CandType const& candidates)
  {
    const auto vtxZ = col.posZ();
    const auto sizeCand = candidates.size();
    if (sizeCand == 0)
      return;

    const auto spher = colCuts.computeSphericity(col, tracks);
    float mult = 0;
    int multNtr = 0;
    if (isRun3) {
      if (useCent) {
        mult = col.centFT0M();
      } else {
        mult = 0;
      }
      multNtr = col.multNTracksPV();
    } else {
      mult = 1; // multiplicity percentile is know in Run 2
      multNtr = col.multTracklets();
    }

    colCuts.fillQA(col, mult);

    // check whether the basic event selection criteria are fulfilled
    // that included checking if there is at least on usable track or V0
    if (!colCuts.isSelectedCollision(col)) {
      return;
    }

    if (colCuts.isEmptyCollision(col, tracks, trackCuts)) {
      return;
    }
    // Filling candidate properties
    rowCandCharmHad.reserve(sizeCand);
    bool isTrackFilled = false;
    for (const auto& candidate : candidates) {
      std::array<float, 3> outputMlPKPi{-1., -1., -1.};
      std::array<float, 3> outputMlPiKP{-1., -1., -1.};
      if constexpr (useCharmMl) {
        /// fill with ML information
        /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
        if (candidate.mlProbLcToPKPi().size() > 0) {
          outputMlPKPi.at(0) = candidate.mlProbLcToPKPi()[0]; /// bkg score
          outputMlPKPi.at(1) = candidate.mlProbLcToPKPi()[1]; /// prompt score
          outputMlPKPi.at(2) = candidate.mlProbLcToPKPi()[2]; /// non-prompt score
        }
        if (candidate.mlProbLcToPiKP().size() > 0) {
          outputMlPiKP.at(0) = candidate.mlProbLcToPiKP()[0]; /// bkg score
          outputMlPiKP.at(1) = candidate.mlProbLcToPiKP()[1]; /// prompt score
          outputMlPiKP.at(2) = candidate.mlProbLcToPiKP()[2]; /// non-prompt score
        }
      }
      auto trackPos1 = candidate.template prong0_as<FemtoHFTracks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<FemtoHFTracks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<FemtoHFTracks>(); // positive daughter (negative for the antiparticles)

      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float BDTScoreBkg,
                           float BDTScorePrompt,
                           float BDTScoreFD) {
        if (FunctionSelection >= 1){
        // Fill tracks if it is not filled for Lc Candidate in an event
            if (!isTrackFilled) {
                isTrackFilled = fillTracksForCharmHadron<false>(tracks, trackPos1, trackNeg, trackPos2, sizeCand);

                // If track filling was successful, fill the collision table
                if (isTrackFilled) {
                    outputCollision(vtxZ, mult, multNtr, spher, magField);
                }
            }

            // fill collision table if track table is filled, i.e., there is at least one Lc-p pair
            if (isTrackFilled) {
                // Row for candidate charm hadron
                rowCandCharmHad(
                    outputCollision.lastIndex(),
                    trackPos1.sign() + trackNeg.sign() + trackPos2.sign(),
                    trackPos1.globalIndex(),
                    trackNeg.globalIndex(),
                    trackPos2.globalIndex(),
                    trackPos1.pt(),
                    trackNeg.pt(),
                    trackPos2.pt(),
                    trackPos1.eta(),
                    trackNeg.eta(),
                    trackPos2.eta(),
                    trackPos1.phi(),
                    trackNeg.phi(),
                    trackPos2.phi(),
                    1 << CandFlag,
                    BDTScoreBkg,
                    BDTScorePrompt,
                    BDTScoreFD);

                // Row for MC candidate charm hadron (if constexpr isMC)
                if constexpr (isMC) {
                    rowCandMCCharmHad(
                        candidate.flagMcMatchRec(),
                        candidate.originMcRec());
                }
            }
      } };

      fillTable(0, candidate.isSelLcToPKPi(), outputMlPKPi.at(0), outputMlPKPi.at(1), outputMlPKPi.at(2));
      fillTable(1, candidate.isSelLcToPiKP(), outputMlPiKP.at(0), outputMlPiKP.at(1), outputMlPiKP.at(2));
    }
  }

  template <typename ParticleType>
  void fillCharmHadMcGen(ParticleType particles)
  {
    // Filling particle properties
    rowCandCharmHadGen.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        rowCandCharmHadGen(
          particle.mcCollisionId(),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  void processDataCharmHad(FemtoFullCollision const& col,
                           aod::BCsWithTimestamps const&,
                           FemtoHFTracks const& tracks,
                           soa::Filtered<CandidateLc> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducer, processDataCharmHad,
                 "Provide experimental data for charm hadron femto", false);

  void processDataCharmHadWithML(FemtoFullCollision const& col,
                                 aod::BCsWithTimestamps const&,
                                 FemtoHFTracks const& tracks,
                                 soa::Filtered<soa::Join<CandidateLc,
                                                         aod::HfMlLcToPKPi>> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducer, processDataCharmHadWithML,
                 "Provide experimental data for charm hadron femto with ml", false);

  void processMCCharmHad(FemtoFullCollisionMC const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<FemtoHFTracks,
                                   aod::McTrackLabels> const& tracks,
                         soa::Filtered<CandidateLcMC> const& candidates,
                         GeneratedMC const& /*particles*/)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducer, processMCCharmHad, "Provide MC for charm hadron", false);

  void processMCCharmHadWithML(FemtoFullCollisionMC const& col,
                               aod::BCsWithTimestamps const&,
                               soa::Join<FemtoHFTracks,
                                         aod::McTrackLabels> const& tracks,
                               soa::Filtered<soa::Join<CandidateLcMC,
                                                       aod::HfMlLcToPKPi>> const& candidates,
                               GeneratedMC const& /*particles*/)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoDreamProducer, processMCCharmHadWithML, "Provide MC for charm hadron with ml", false);

  void processMCCharmHadGen(GeneratedMC const& particles)
  {

    fillCharmHadMcGen(particles);
  }
  PROCESS_SWITCH(femtoDreamProducer, processMCCharmHadGen, "Provide MC Generated charm hadron", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<femtoDreamProducer>(cfgc)};
}
