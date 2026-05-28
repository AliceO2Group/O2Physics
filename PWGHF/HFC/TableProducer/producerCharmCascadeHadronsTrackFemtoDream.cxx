// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file producerCharmCascadeHadronsTrackFemtoDream.cxx
/// \brief FemtoDream producer for Ξc± → Ξ∓ π± π± + associated track (e.g. π)
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfMlResponseXicToXiPiPi.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;
using namespace o2::aod::hf_cand_xic_to_xi_pi_pi;

enum Event : uint8_t {
  All = 0,
  RejEveSel,
  RejNoTracksAndCharm,
  TrackSelected,
  CharmSelected,
  PairSelected
};

enum MlMode : uint8_t {
  NoMl = 0,
  FillMlFromSelector,
  FillMlFromNewBDT
};

struct HfProducerCharmCascadeHadronsTrackFemtoDream {

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDColMasks> rowMasks;
  Produces<aod::FDHfCand3Prong> rowCandCharm3Prong;
  Produces<aod::FDHfCand3ProngXic> rowCandCharm3ProngXic;
  Produces<aod::FDHfCandMC> rowCandMcCharmHad;
  Produces<aod::FDHfCandMCGen> rowCandCharmHadGen;
  Produces<aod::FDParticlesIndex> outputPartsIndex;
  Produces<aod::FDTrkTimeStamp> outputPartsTime;
  Produces<aod::FDMCCollisions> outputMcCollision;
  Produces<aod::FDMCCollLabels> outputCollsMcLabels;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMc;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMcLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMc;
  Produces<aod::FDExtMCLabels> outputPartsExtMcLabels;

  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTXicToXiPiPi"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_XicToXiPiPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  } ccdbCfg;

  struct : ConfigurableGroup {
    Configurable<int> applyMlMode{"applyMlMode", 1, "None: 0, BDT model from candidate selector: 1, New BDT model on Top of candidate selector: 2"};
    Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
    Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
    Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
    Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  } mlCfg;

  struct : ConfigurableGroup {
    Configurable<float> pTrackMethod1Max{"pTrackMethod1Max", 0.85f, "Kaon PID Method1 (TPC-only): maximum p (GeV/c)"};
    Configurable<float> pTrackExcludeMin{"pTrackExcludeMin", 0.50f, "Kaon PID Method1: excluded p window minimum (GeV/c)"};
    Configurable<float> pTrackExcludeMax{"pTrackExcludeMax", 0.65f, "Kaon PID Method1: excluded p window maximum (GeV/c)"};
    Configurable<float> pTrackPiRejMin{"pTrackPiRejMin", 0.50f, "Kaon PID Method1: pion rejection active for p > this (GeV/c)"};
    Configurable<float> pTrackElRejMin{"pTrackElRejMin", 0.30f, "Kaon PID Method1: electron rejection active for p > this (GeV/c)"};
    Configurable<float> pTrackTightMin{"pTrackTightMin", 1.20f, "Kaon PID Method2 (TPC+TOF): tighten cuts for p > this (GeV/c)"};
    Configurable<float> nSigmaTpcKaMax{"nSigmaTpcKaMax", 3.f, "Kaon PID Method1: require |nSigmaTpcKa| < this"};
    Configurable<float> nSigmaTpcPiMin{"nSigmaTpcPiMin", 3.f, "Kaon PID Method1: require |nSigmaTpcPi| > this (pion)"};
    Configurable<float> nSigmaTpcElMin{"nSigmaTpcElMin", 3.f, "Kaon PID Method1: require |nSigmaTpcEl| > this (electron)"};
    Configurable<float> nSigmaCombKaMax{"nSigmaCombKaMax", 3.f, "Kaon PID Method2: require |nSigmaCombKa| < this"};
    Configurable<float> nSigmaCombKaTightMax{"nSigmaCombKaTightMax", 2.f, "Kaon PID Method2: for p > pTrackTightMin require |nSigmaCombKa| < this"};
    Configurable<float> nSigmaCombPiMax{"nSigmaCombPiMax", 6.f, "Kaon PID Method2: for p > pTrackTightMin require |nSigmaCombPi| < this"};
  } kaonPidSel;

  Configurable<bool> isDebug{"isDebug", true, "Enable Debug tables"};
  Configurable<bool> isRun3{"isRun3", true, "Running on Run3 or pilot"};
  Configurable<int> selectionFlagHadron{"selectionFlagHadron", 1, "Selection flag for Ξc± → Ξππ (same bit mask as hf_sel_candidate_xic::isSelXicToXiPiPi)"};
  Configurable<bool> useCent{"useCent", false, "Enable centrality for charm hadron"};

  Configurable<int> trkPDGCode{"trkPDGCode", kPiPlus, "PDG code of the associated track for MC truth and PID (default: π)"};
  Configurable<std::vector<float>> trkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "trk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "trk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "trk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> trkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "trk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
  Configurable<std::vector<float>> trkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "trk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> trkPIDnSigmaOffsetTPC{"trkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> trkPIDnSigmaOffsetTOF{"trkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<float>> trkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "trk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> trkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "trk"), std::vector<float>{0.15f, 0.4f, 0.6f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "trk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "trk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "trk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> trkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "trk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "trk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> trkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "trk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};

  FemtoDreamTrackSelection trackCuts;
  o2::analysis::HfMlResponseXicToXiPiPi<float> hfMlResponseXic;
  std::vector<float> outputMlXic;
  o2::ccdb::CcdbApi ccdbApi;
  o2::hf_evsel::HfEventSelection hfEvSel;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};

  float magField{};
  int runNumber{};

  using CandidateXic = soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>;
  using CandidateXicMc = soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfCandXicMcRec>;
  using CandidateXicKf = soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfSelXicToXiPiPi>;
  using CandidateXicKfMc = soa::Join<aod::HfCandXic, aod::HfCandXicKF, aod::HfSelXicToXiPiPi, aod::HfCandXicMcRec>;

  using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
  using FemtoFullCollisionMc = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
  using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using FemtoHFMcTracks = soa::Join<aod::McTrackLabels, FemtoHFTracks>;
  using GeneratedXicMc = soa::Join<aod::McParticles, aod::HfCandXicMcGen>;

  Filter filterSelectCandidateXic = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagHadron;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry trackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext&)
  {
    std::array<bool, 9> processes = {doprocessDataXicToXiPiPi, doprocessDataXicToXiPiPiKf,
                                     doprocessDataXicToXiPiPiWithML, doprocessDataXicToXiPiPiWithMLKf,
                                     doprocessMcXicToXiPiPi, doprocessMcXicToXiPiPiKf,
                                     doprocessMcXicToXiPiPiWithML, doprocessMcXicToXiPiPiWithMLKf,
                                     doprocessMcXicToXiPiPiGen};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int const cutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    trackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});

    constexpr int kEventTypes = PairSelected + 1;
    std::string labels[kEventTypes];
    labels[Event::All] = "All events";
    labels[Event::RejEveSel] = "rejected by event selection";
    labels[Event::RejNoTracksAndCharm] = "rejected by no tracks and charm";
    labels[Event::TrackSelected] = "with tracks ";
    labels[Event::CharmSelected] = "with charm hadrons ";
    labels[Event::PairSelected] = "with pairs";

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
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &trackRegistry);

    runNumber = 0;
    magField = 0.f;
    ccdb->setURL(ccdbCfg.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    hfEvSel.init(qaRegistry, &zorroSummary);

    int64_t const now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    bool useXicMl = doprocessDataXicToXiPiPiWithML || doprocessMcXicToXiPiPiWithML || doprocessDataXicToXiPiPiWithMLKf || doprocessMcXicToXiPiPiWithMLKf;
    if (mlCfg.applyMlMode == FillMlFromNewBDT && useXicMl) {
      hfMlResponseXic.configure(mlCfg.binsPtMl, mlCfg.cutsMl, mlCfg.cutDirMl, mlCfg.nClassesMl);
      hfMlResponseXic.cacheInputFeaturesIndices(mlCfg.namesInputFeatures);
      if (ccdbCfg.loadModelsFromCCDB) {
        ccdbApi.init(ccdbCfg.ccdbUrl);
        hfMlResponseXic.setModelPathsCCDB(ccdbCfg.onnxFileNames, ccdbApi, ccdbCfg.modelPathsCCDB, ccdbCfg.timestampCCDB);
      } else {
        hfMlResponseXic.setModelPathsLocal(ccdbCfg.onnxFileNames);
      }
      hfMlResponseXic.init();
    }
  }

  void getMagneticFieldTesla(const aod::BCsWithTimestamps::iterator& bc)
  {
    initCCDB(bc, runNumber, ccdb, !isRun3 ? ccdbCfg.ccdbPathGrp : ccdbCfg.ccdbPathGrpMag, lut, !isRun3);
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
                     -999., -999., -999., -999.,
                     -999., -999., -999., -999.,
                     -999., -999., -999.);
  }

  template <typename CollisionType, typename ParticleType>
  void fillMcParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      auto particleMc = particle.mcParticle();
      auto pdgCode = particleMc.pdgCode();
      int particleOrigin = 99;
      int pdgCodeMother = -1;
      constexpr int GenFromTransport = -1;
      auto motherparticlesMc = particleMc.template mothers_as<aod::McParticles>();
      if (std::abs(pdgCode) == std::abs(trkPDGCode.value)) {
        if ((col.has_mcCollision() && (particleMc.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
        } else if (particleMc.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (particleMc.getProcess() == TMCProcess::kPDecay && particleMc.getGenStatusCode() == GenFromTransport && !motherparticlesMc.empty()) {
          auto motherparticleMc = motherparticlesMc.front();
          pdgCodeMother = motherparticleMc.pdgCode();
          particleOrigin = checkDaughterType(fdparttype, motherparticleMc.pdgCode());
        } else if (particleMc.getProcess() == TMCProcess::kPHInhelastic && particleMc.getGenStatusCode() == GenFromTransport) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
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

  template <typename TrackType>
  bool isTrackKaonPidSelected(const TrackType& track)
  {
    const float pTrack = track.p();

    bool isTrackKaonPidMethod1 = true;
    if (pTrack >= kaonPidSel.pTrackMethod1Max) {
      isTrackKaonPidMethod1 = false;
    }
    if (std::abs(track.tpcNSigmaKa()) >= kaonPidSel.nSigmaTpcKaMax) {
      isTrackKaonPidMethod1 = false;
    }
    if (pTrack >= kaonPidSel.pTrackExcludeMin && pTrack <= kaonPidSel.pTrackExcludeMax) {
      isTrackKaonPidMethod1 = false;
    }
    if (pTrack > kaonPidSel.pTrackPiRejMin && std::abs(track.tpcNSigmaPi()) <= kaonPidSel.nSigmaTpcPiMin) {
      isTrackKaonPidMethod1 = false;
    }
    if (pTrack > kaonPidSel.pTrackElRejMin && std::abs(track.tpcNSigmaEl()) <= kaonPidSel.nSigmaTpcElMin) {
      isTrackKaonPidMethod1 = false;
    }

    bool isTrackKaonPidMethod2 = true;
    if (pTrack > kaonPidSel.pTrackMethod1Max && !track.hasTOF()) {
      isTrackKaonPidMethod2 = false;
    }

    const float nSigmaCombKa = std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa());
    const float nSigmaCombPi = std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi());

    if (std::abs(nSigmaCombKa) >= kaonPidSel.nSigmaCombKaMax) {
      isTrackKaonPidMethod2 = false;
    }

    if (pTrack > kaonPidSel.pTrackTightMin) {
      if (std::abs(nSigmaCombKa) >= kaonPidSel.nSigmaCombKaTightMax) {
        isTrackKaonPidMethod2 = false;
      }
      if (std::abs(nSigmaCombPi) <= kaonPidSel.nSigmaCombPiMax) {
        isTrackKaonPidMethod2 = false;
      }
    }

    return isTrackKaonPidMethod1 || isTrackKaonPidMethod2;
  }

  template <typename CollisionType>
  void fillMcCollision(CollisionType const& col)
  {
    if (col.has_mcCollision()) {
      outputCollsMcLabels(outputMcCollision.lastIndex());
    } else {
      outputCollsMcLabels(-1);
    }
  }

  template <bool IsMc = false, typename TrackType, typename CollisionType>
  bool fillTracksForCharmHadron(CollisionType const& col, TrackType const& tracks)
  {
    std::vector<int> const childIDs = {0, 0};
    bool fIsTrackFilled = false;

    for (const auto& track : tracks) {
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, true>(track);
      auto cutContainer = trackCuts.getCutContainer<false, aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));
      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      int64_t timeStamp = bc.timestamp();
      outputPartsIndex(track.globalIndex());
      outputPartsTime(timeStamp);

      if (trkPDGCode == kKPlus) {
        const auto pidTrackPassBit = static_cast<aod::femtodreamparticle::cutContainerType>(isTrackKaonPidSelected(track));

        outputParts(outputCollision.lastIndex(),
                    track.pt(),
                    track.eta(),
                    track.phi(),
                    aod::femtodreamparticle::ParticleType::kTrack,
                    cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                    pidTrackPassBit,
                    track.dcaXY(), childIDs, 0, 0);
      } else {
        outputParts(outputCollision.lastIndex(),
                    track.pt(),
                    track.eta(),
                    track.phi(),
                    aod::femtodreamparticle::ParticleType::kTrack,
                    cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                    cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                    track.dcaXY(), childIDs, 0, 0);
      }

      fIsTrackFilled = true;
      if (isDebug.value) {
        fillDebugParticle(track);
      }

      if constexpr (IsMc) {
        fillMcParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }
    return fIsTrackFilled;
  }

  template <typename C, typename T, typename TC>
  bool isNoSelectedTracks(C const& /*col*/, T const& tracks, TC& cuts)
  {
    for (auto const& track : tracks) {
      if (cuts.isSelectedMinimal(track)) {
        return false;
      }
    }
    return true;
  }

  template <bool IsMc, bool UseCharmMl, typename TrackType, typename CollisionType, typename CandType>
  void fillXicHadronTable(CollisionType const& col, TrackType const& tracks, CandType const& candidates)
  {
    const auto vtxZ = col.posZ();
    const auto sizeCand = candidates.size();
    const auto spher = 2.f;
    float mult = 0;
    int multNtr = 0;
    if (isRun3) {
      mult = useCent ? col.centFT0M() : 0.f;
      multNtr = col.multNTracksPV();
    } else {
      mult = 1.f;
      multNtr = col.multTracklets();
    }

    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(col, mult, ccdb, qaRegistry);
    qaRegistry.fill(HIST("hEventQA"), 1 + Event::All);
    hfEvSel.fillHistograms(col, rejectionMask, mult);
    if (rejectionMask != 0) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::RejEveSel);
      return;
    }

    if (isNoSelectedTracks(col, tracks, trackCuts) && sizeCand <= 0) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::RejNoTracksAndCharm);
      return;
    }

    outputCollision(vtxZ, mult, multNtr, spher, magField);
    if constexpr (IsMc) {
      fillMcCollision(col);
    }

    rowCandCharm3Prong.reserve(sizeCand);
    rowCandCharm3ProngXic.reserve(sizeCand);
    bool isTrackFilled = false;
    int nSelectedXic = 0;

    for (const auto& candidate : candidates) {
      if (candidate.isSelXicToXiPiPi() < selectionFlagHadron) {
        continue;
      }

      outputMlXic = {-1.f, -1.f, -1.f};
      bool isSelectedMlXicToXiPiPi = true;
      if constexpr (UseCharmMl) {
        if (mlCfg.applyMlMode == FillMlFromSelector) {
          if (candidate.mlProbXicToXiPiPi().size() > 0) {
            outputMlXic.at(0) = candidate.mlProbXicToXiPiPi()[0];
            outputMlXic.at(1) = candidate.mlProbXicToXiPiPi()[1];
            outputMlXic.at(2) = candidate.mlProbXicToXiPiPi()[2];
          }
        } else if (mlCfg.applyMlMode == FillMlFromNewBDT) {
          isSelectedMlXicToXiPiPi = false;
          if (candidate.mlProbXicToXiPiPi().size() > 0) {
            std::vector<float> inputFeaturesXicToXiPiPi = hfMlResponseXic.getInputFeatures(candidate);
            isSelectedMlXicToXiPiPi = hfMlResponseXic.isSelectedMl(inputFeaturesXicToXiPiPi, candidate.pt(), outputMlXic);
          }
          if (!isSelectedMlXicToXiPiPi) {
            continue;
          }
        } else {
          LOGF(fatal, "Please check your ML configuration.");
        }
      }

      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      int64_t timeStamp = bc.timestamp();
      const auto eta0 = static_cast<float>(RecoDecay::eta(std::array{candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()}));
      const auto eta1 = static_cast<float>(RecoDecay::eta(std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()}));
      const auto eta2 = static_cast<float>(RecoDecay::eta(std::array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()}));
      const auto phi0 = static_cast<float>(RecoDecay::phi(candidate.pxProng0(), candidate.pyProng0()));
      const auto phi1 = static_cast<float>(RecoDecay::phi(candidate.pxProng1(), candidate.pyProng1()));
      const auto phi2 = static_cast<float>(RecoDecay::phi(candidate.pxProng2(), candidate.pyProng2()));

      rowCandCharm3Prong(
        outputCollision.lastIndex(),
        timeStamp,
        candidate.sign(),
        candidate.pi0Id(),
        candidate.pi1Id(),
        candidate.bachelorId(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        eta0,
        eta1,
        eta2,
        phi0,
        phi1,
        phi2,
        1 << 0,
        outputMlXic.at(0),
        outputMlXic.at(1),
        outputMlXic.at(2));

      rowCandCharm3ProngXic(
        candidate.invMassXicPlus(),
        candidate.posTrackId(),
        candidate.negTrackId());

      ++nSelectedXic;
      if constexpr (IsMc) {
        rowCandMcCharmHad(
          candidate.flagMcMatchRec(),
          candidate.originMcRec());
      }
    }

    isTrackFilled = fillTracksForCharmHadron<IsMc>(col, tracks);

    aod::femtodreamcollision::BitMaskType bitTrack = 0;
    if (isTrackFilled) {
      bitTrack |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::TrackSelected);
    }

    aod::femtodreamcollision::BitMaskType bitCand = 0;
    if (nSelectedXic > 0) {
      bitCand |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::CharmSelected);
    }

    if (isTrackFilled && nSelectedXic > 0) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::PairSelected);
    }

    rowMasks(bitTrack, bitCand, 0);
  }

  void fillXicMcGen(GeneratedXicMc const& particles)
  {
    rowCandCharmHadGen.reserve(particles.size());
    for (const auto& particle : particles) {
      const int absFlag = std::abs(static_cast<int>(particle.flagMcMatchGen()));
      if (absFlag == (1 << DecayType::XicToXiPiPi) || absFlag == (1 << DecayType::XicToXiResPiToXiPiPi)) {
        rowCandCharmHadGen(
          particle.mcCollisionId(),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  void processDataXicToXiPiPi(FemtoFullCollision const& col,
                              aod::BCsWithTimestamps const&,
                              FemtoHFTracks const& tracks,
                              soa::Filtered<CandidateXic> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processDataXicToXiPiPi, "Data for XicToXiPiPi femto (DCAFitter; no HFCANDXICKF)", false);

  void processDataXicToXiPiPiKf(FemtoFullCollision const& col,
                                aod::BCsWithTimestamps const&,
                                FemtoHFTracks const& tracks,
                                soa::Filtered<CandidateXicKf> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processDataXicToXiPiPiKf, "Data for XicToXiPiPi femto (KFParticle; requires HFCANDXICKF)", false);

  void processDataXicToXiPiPiWithML(FemtoFullCollision const& col,
                                    aod::BCsWithTimestamps const&,
                                    FemtoHFTracks const& tracks,
                                    soa::Filtered<soa::Join<CandidateXic, aod::HfMlXicToXiPiPi>> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processDataXicToXiPiPiWithML, "Data for XicToXiPiPi femto with ML (DCAFitter)", false);

  void processDataXicToXiPiPiWithMLKf(FemtoFullCollision const& col,
                                      aod::BCsWithTimestamps const&,
                                      FemtoHFTracks const& tracks,
                                      soa::Filtered<soa::Join<CandidateXicKf, aod::HfMlXicToXiPiPi>> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processDataXicToXiPiPiWithMLKf, "Data for XicToXiPiPi femto with ML (KFParticle)", false);

  void processMcXicToXiPiPi(FemtoFullCollisionMc const& col,
                          aod::BCsWithTimestamps const&,
                          FemtoHFMcTracks const& tracks,
                          aod::McParticles const&,
                          soa::Filtered<CandidateXicMc> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processMcXicToXiPiPi, "MC for XicToXiPiPi (DCAFitter)", false);

  void processMcXicToXiPiPiKf(FemtoFullCollisionMc const& col,
                             aod::BCsWithTimestamps const&,
                             FemtoHFMcTracks const& tracks,
                             aod::McParticles const&,
                             soa::Filtered<CandidateXicKfMc> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processMcXicToXiPiPiKf, "MC for XicToXiPiPi (KFParticle)", false);

  void processMcXicToXiPiPiWithML(FemtoFullCollisionMc const& col,
                                aod::BCsWithTimestamps const&,
                                FemtoHFMcTracks const& tracks,
                                aod::McParticles const&,
                                soa::Filtered<soa::Join<CandidateXicMc, aod::HfMlXicToXiPiPi>> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processMcXicToXiPiPiWithML, "MC for XicToXiPiPi with ML (DCAFitter)", false);

  void processMcXicToXiPiPiWithMLKf(FemtoFullCollisionMc const& col,
                                  aod::BCsWithTimestamps const&,
                                  FemtoHFMcTracks const& tracks,
                                  aod::McParticles const&,
                                  soa::Filtered<soa::Join<CandidateXicKfMc, aod::HfMlXicToXiPiPi>> const& candidates)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillXicHadronTable<true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processMcXicToXiPiPiWithMLKf, "MC for XicToXiPiPi with ML (KFParticle)", false);

  void processMcXicToXiPiPiGen(GeneratedXicMc const& particles)
  {
    fillXicMcGen(particles);
  }
  PROCESS_SWITCH(HfProducerCharmCascadeHadronsTrackFemtoDream, processMcXicToXiPiPiGen, "MC generated XicToXiPiPi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfProducerCharmCascadeHadronsTrackFemtoDream>(cfgc)};
}
