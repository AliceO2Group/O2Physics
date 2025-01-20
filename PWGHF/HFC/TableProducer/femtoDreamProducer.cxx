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
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include <string>
#include <vector>
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
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Core/CentralityEstimation.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;

// event types
enum Event : uint8_t {
  kAll = 0,
  kRejEveSel,
  kRejNoTracksAndCharm,
  kTrackSelected,
  kCharmSelected,
  kPairSelected
};

struct HfFemtoDreamProducer {

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDColMasks> rowMasks;
  Produces<aod::FDHfCand> rowCandCharmHad;
  Produces<aod::FDHfCandMC> rowCandMcCharmHad;
  Produces<aod::FDHfCandMCGen> rowCandCharmHadGen;
  Produces<aod::FDParticlesIndex> outputPartsIndex;
  Produces<aod::FDMCCollisions> outputMcCollision;
  Produces<aod::FDMCCollLabels> outputCollsMcLabels;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMc;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMcLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMc;
  Produces<aod::FDExtMCLabels> outputPartsExtMcLabels;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // Configurable<bool> isForceGRP{"isForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  Configurable<bool> isDebug{"isDebug", true, "Enable Debug tables"};
  Configurable<bool> isRun3{"isRun3", true, "Running on Run3 or pilot"};

  /// Lc table
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<bool> useCent{"useCent", false, "Enable centrality for lc"};

  Configurable<int> trkPDGCode{"trkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  Configurable<std::vector<float>> trkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "trk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "trk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "trk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> trkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "trk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton}, "Trk sel: Particles species for PID"};
  Configurable<std::vector<float>> trkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "trk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> trkPIDnSigmaOffsetTPC{"trkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> trkPIDnSigmaOffsetTOF{"trkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<float>> trkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "trk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> trkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "trk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "trk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "trk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "trk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "trk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "trk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "trk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};

  using CandidateLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;
  using CandidateLcMc = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>;

  using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
  using FemtoFullCollisionMc = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
  using FemtoFullMcgenCollisions = soa::Join<aod::McCollisions, o2::aod::MultsExtraMC>;
  using FemtoFullMcgenCollision = FemtoFullMcgenCollisions::iterator;
  using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
  using FemtoHFTrack = FemtoHFTracks::iterator;
  using FemtoHFMcTracks = soa::Join<aod::McTrackLabels, FemtoHFTracks>;
  using FemtoHFMcTrack = FemtoHFMcTracks::iterator;

  using GeneratedMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;

  FemtoDreamTrackSelection trackCuts;

  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry TrackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};

  HfHelper hfHelper;
  o2::hf_evsel::HfEventSelection hfEvSel;

  float magField;
  int runNumber;

  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  o2::base::MatLayerCylSet* lut;
  // if (doPvRefit){ lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));} //! may be it useful, will check later

  void init(InitContext&)
  {
    std::array<bool, 5> processes = {doprocessDataCharmHad, doprocessMcCharmHad, doprocessDataCharmHadWithML, doprocessMcCharmHadWithML, doprocessMcCharmHadGen};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int CutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    TrackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{CutBits + 1, -0.5, CutBits + 0.5}});

    // event QA histograms
    constexpr int kEventTypes = kPairSelected + 1;
    std::string labels[kEventTypes];
    labels[Event::kAll] = "All events";
    labels[Event::kRejEveSel] = "rejected by event selection";
    labels[Event::kRejNoTracksAndCharm] = "rejected by no tracks and charm";
    labels[Event::kTrackSelected] = "with tracks ";
    labels[Event::kCharmSelected] = "with charm hadrons ";
    labels[Event::kPairSelected] = "with pairs";

    static const AxisSpec axisEvents = {kEventTypes, 0.5, kEventTypes + 0.5, ""};
    qaRegistry.add("hEventQA", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kEventTypes; iBin++) {
      qaRegistry.get<TH1>(HIST("hEventQA"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

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

    hfEvSel.addHistograms(qaRegistry); // collision monitoring

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    initCCDB(bc, runNumber, ccdb, !isRun3 ? ccdbPathGrp : ccdbPathGrpMag, lut, !isRun3);
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
                     -999.,
                     particle.tpcNSigmaPi(),
                     particle.tpcNSigmaKa(),
                     particle.tpcNSigmaPr(),
                     particle.tpcNSigmaDe(),
                     -999.,
                     -999.,
                     -999.,
                     particle.tofNSigmaPi(),
                     particle.tofNSigmaKa(),
                     particle.tofNSigmaPr(),
                     particle.tofNSigmaDe(),
                     -999., -999., -999., -999.,
                     -999., -999., -999., -999.,
                     -999., -999., -999., -999.,
                     -999., -999., -999., -999.);
  }

  template <typename CollisionType, typename ParticleType>
  void fillMcParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMc = particle.mcParticle();
      auto pdgCode = particleMc.pdgCode();
      int particleOrigin = 99;
      int pdgCodeMother = -1;
      // get list of mothers, but it could be empty (for example in case of injected light nuclei)
      auto motherparticlesMc = particleMc.template mothers_as<aod::McParticles>();
      // check pdg code
      // if this fails, the particle is a fake
      if (abs(pdgCode) == abs(trkPDGCode.value)) {
        // check first if particle is from pile up
        // check if the collision associated with the particle is the same as the analyzed collision by checking their Ids
        if ((col.has_mcCollision() && (particleMc.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
          // check if particle is primary
        } else if (particleMc.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
          // check if particle is secondary
          // particle is from a decay -> getProcess() == 4
          // particle is generated during transport -> getGenStatusCode() == -1
          // list of mothers is not empty
        } else if (particleMc.getProcess() == 4 && particleMc.getGenStatusCode() == -1 && !motherparticlesMc.empty()) {
          // get direct mother
          auto motherparticleMc = motherparticlesMc.front();
          pdgCodeMother = motherparticleMc.pdgCode();
          particleOrigin = checkDaughterType(fdparttype, motherparticleMc.pdgCode());
          // check if particle is material
          // particle is from inelastic hadronic interaction -> getProcess() == 23
          // particle is generated during transport -> getGenStatusCode() == -1
        } else if (particleMc.getProcess() == 23 && particleMc.getGenStatusCode() == -1) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
          // cross check to see if we missed a case
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
        // if pdg code is wrong, particle is fake
      } else {
        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMc(particleOrigin, pdgCode, particleMc.pt(), particleMc.eta(), particleMc.phi());
      outputPartsMcLabels(outputPartsMc.lastIndex());
      if (isDebug) {
        outputPartsExtMcLabels(outputPartsMc.lastIndex());
        outputDebugPartsMc(pdgCodeMother);
      }
    } else {
      outputPartsMcLabels(-1);
      if (isDebug) {
        outputPartsExtMcLabels(-1);
      }
    }
  }

  template <typename CollisionType>
  void fillMcCollision(CollisionType const& col)
  {
    if (col.has_mcCollision()) {
      // auto genMCcol = col.template mcCollision_as<FemtoFullMcgenCollisions>();
      // outputMcCollision(genMCcol.multMCNParticlesEta08());
      outputCollsMcLabels(outputMcCollision.lastIndex());
    } else {
      outputCollsMcLabels(-1);
    }
  }

  template <bool isMc = false, typename TrackType, typename CollisionType>
  bool fillTracksForCharmHadron(CollisionType const& col, TrackType const& tracks)
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

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, true>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<false, aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

      // track global index
      outputPartsIndex(track.globalIndex());
      // now the table is filled

      outputParts(outputCollision.lastIndex(),
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

      if constexpr (isMc) {
        fillMcParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }
    return fIsTrackFilled;
  }

  template <bool isMc, bool useCharmMl, typename TrackType, typename CollisionType, typename CandType>
  void fillCharmHadronTable(CollisionType const& col, TrackType const& tracks, CandType const& candidates)
  {
    const auto vtxZ = col.posZ();
    const auto sizeCand = candidates.size();
    const auto spher = 2.; // dummy value for the moment
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

    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(col, mult, ccdb, qaRegistry);

    qaRegistry.fill(HIST("hEventQA"), 1 + Event::kAll);

    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(col, rejectionMask, mult);
    if (rejectionMask != 0) {
      /// at least one event selection not satisfied --> reject the candidate
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::kRejEveSel);
      return;
    }

    if (isNoSelectedTracks(col, tracks, trackCuts) && sizeCand <= 0) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::kRejNoTracksAndCharm);
      return;
    }

    outputCollision(vtxZ, mult, multNtr, spher, magField);
    if constexpr (isMc) {
      fillMcCollision(col);
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
      auto trackPos1 = candidate.template prong0_as<TrackType>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<TrackType>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<TrackType>(); // positive daughter (negative for the antiparticles)

      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float BDTScoreBkg,
                           float BDTScorePrompt,
                           float BDTScoreFD) {
        if (FunctionSelection >= 1){
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

                // Row for MC candidate charm hadron (if constexpr isMc)
                if constexpr (isMc) {
                  rowCandMcCharmHad(
                    candidate.flagMcMatchRec(),
                    candidate.originMcRec());}
      } };

      fillTable(0, candidate.isSelLcToPKPi(), outputMlPKPi.at(0), outputMlPKPi.at(1), outputMlPKPi.at(2));
      fillTable(1, candidate.isSelLcToPiKP(), outputMlPiKP.at(0), outputMlPiKP.at(1), outputMlPiKP.at(2));
    }

    if (!isTrackFilled) {
      isTrackFilled = fillTracksForCharmHadron<isMc>(col, tracks);
      // If track filling was successful, fill the collision table
    }

    aod::femtodreamcollision::BitMaskType bitTrack = 0;
    if (isTrackFilled) {
      bitTrack |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::kTrackSelected);
    }

    aod::femtodreamcollision::BitMaskType bitCand = 0;
    if (sizeCand > 0) {
      bitCand |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::kCharmSelected);
    }

    if (isTrackFilled && (sizeCand > 0))
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::kPairSelected);

    rowMasks(static_cast<aod::femtodreamcollision::BitMaskType>(bitTrack),
             static_cast<aod::femtodreamcollision::BitMaskType>(bitCand),
             0);
  }

  // check if there is no selected track
  /// \param C type of the collision
  /// \param T type of the tracks
  /// \param TC type of the femto track cuts
  /// \return whether or not the tracks fulfills the all selections
  template <typename C, typename T, typename TC>
  bool isNoSelectedTracks(C const& /*col*/, T const& tracks, TC& trackCuts)
  {
    for (auto const& track : tracks) {
      if (trackCuts.isSelectedMinimal(track)) {
        return false;
      }
    }
    return true;
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
  PROCESS_SWITCH(HfFemtoDreamProducer, processDataCharmHad,
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
  PROCESS_SWITCH(HfFemtoDreamProducer, processDataCharmHadWithML,
                 "Provide experimental data for charm hadron femto with ml", false);

  void processMcCharmHad(FemtoFullCollisionMc const& col,
                         aod::BCsWithTimestamps const&,
                         FemtoHFMcTracks const& tracks,
                         aod::McParticles const&,
                         CandidateLcMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfFemtoDreamProducer, processMcCharmHad, "Provide Mc for charm hadron", false);

  void processMcCharmHadWithML(FemtoFullCollisionMc const& col,
                               aod::BCsWithTimestamps const&,
                               FemtoHFMcTracks const& tracks,
                               aod::McParticles const&,
                               soa::Join<CandidateLcMc,
                                         aod::HfMlLcToPKPi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());

    fillCharmHadronTable<true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfFemtoDreamProducer, processMcCharmHadWithML, "Provide Mc for charm hadron with ml", false);

  void processMcCharmHadGen(GeneratedMc const& particles)
  {

    fillCharmHadMcGen(particles);
  }
  PROCESS_SWITCH(HfFemtoDreamProducer, processMcCharmHadGen, "Provide Mc Generated charm hadron", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfFemtoDreamProducer>(cfgc)};
}
