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

/// \file producerCharmHadronsTrackFemtoDream.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Ravindra Singh, GSI, ravindra.singh@cern.ch
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch
/// \author Yunfan Liu, Central China Normal University, yunfan.l@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfMlResponseD0ToKPi.h"
#include "PWGHF/Core/HfMlResponseDplusToPiKPi.h"
#include "PWGHF/Core/HfMlResponseDstarToD0Pi.h"
#include "PWGHF/Core/HfMlResponseLcToPKPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

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

// event types
enum Event : uint8_t {
  All = 0,
  RejEveSel,
  RejNoTracksAndCharm,
  TrackSelected,
  CharmSelected,
  PairSelected
};

// ml modes
enum MlMode : uint8_t {
  NoMl = 0,
  FillMlFromSelector,
  FillMlFromNewBDT
};

// decay channels
enum DecayChannel { DplusToPiKPi = 0,
                    LcToPKPi,
                    D0ToPiK,
                    DstarToD0Pi
};

enum class D0CandFlag : uint8_t {
  D0 = 0,
  D0Bar = 1,
  Reflected = 2
};

struct HfProducerCharmHadronsTrackFemtoDream {

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDColMasks> rowMasks;
  Produces<aod::FDHfCand3Prong> rowCandCharm3Prong;
  Produces<aod::FDHfCand2Prong> rowCandCharm2Prong;
  Produces<aod::FDHfCandDstar> rowCandCharmDstar;
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

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTLc"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_LcToPKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // Configurable<bool> isForceGRP{"isForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  Configurable<bool> isDebug{"isDebug", true, "Enable Debug tables"};
  Configurable<bool> isRun3{"isRun3", true, "Running on Run3 or pilot"};

  /// Charm hadron table
  Configurable<int> selectionFlagHadron{"selectionFlagHadron", 1, "Selection Flag for Charm Hadron: 1 for Lc, 7 for Dplus (Topologic and PID cuts)"};
  Configurable<bool> useCent{"useCent", false, "Enable centrality for Charm Hadron"};

  Configurable<int> trkPDGCode{"trkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  Configurable<std::vector<float>> trkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "trk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "trk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> trkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "trk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> trkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "trk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
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
  // ML inference
  Configurable<int> applyMlMode{"applyMlMode", 1, "None: 0, BDT model from candidate selector: 1, New BDT model on Top of candidate selector: 2"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  FemtoDreamTrackSelection trackCuts;

  o2::analysis::HfMlResponseLcToPKPi<float> hfMlResponseLc;
  o2::analysis::HfMlResponseDplusToPiKPi<float> hfMlResponseDplus;
  o2::analysis::HfMlResponseD0ToKPi<float> hfMlResponseD0;
  o2::analysis::HfMlResponseDstarToD0Pi<float> hfMlResponseDstar;

  std::vector<float> outputMlD0;
  std::vector<float> outputMlD0bar;
  std::vector<float> outputMlDstar;
  std::vector<float> outputMlDplus;
  std::vector<float> outputMlPKPi;
  std::vector<float> outputMlPiKP;
  o2::ccdb::CcdbApi ccdbApi;
  o2::hf_evsel::HfEventSelection hfEvSel;
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; /// Accessing the CCDB
  o2::base::MatLayerCylSet* lut{};
  // if (doPvRefit){ lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));} //! may be it useful, will check later

  float magField{};
  int runNumber{};
  using CandidateD0 = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>;
  using CandidateD0Mc = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using CandidateDstar = soa::Join<aod::HfCandDstarsWPid, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi>;
  using CandidateDstarMc = soa::Join<aod::HfCandDstarsWPid, aod::HfD0FromDstar, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>;
  using CandidateDplus = soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi>;
  using CandidateDplusMc = soa::Join<aod::HfCand3ProngWPidPiKa, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>;
  using CandidateLc = soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc>;
  using CandidateLcMc = soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc, aod::HfCand3ProngMcRec>;

  using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
  using FemtoFullCollisionMc = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
  using FemtoFullMcgenCollisions = soa::Join<aod::McCollisions, o2::aod::MultsExtraMC>;
  using FemtoFullMcgenCollision = FemtoFullMcgenCollisions::iterator;
  using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using FemtoHFTrack = FemtoHFTracks::iterator;
  using FemtoHFMcTracks = soa::Join<aod::McTrackLabels, FemtoHFTracks>;
  using FemtoHFMcTrack = FemtoHFMcTracks::iterator;

  using Generated3ProngMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using Generated2ProngMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
  using GeneratedDstarMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandDstarMcGen>>;

  Filter filterSelectCandidateD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagHadron || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagHadron);
  Filter filterSelectCandidateDstar = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Filter filterSelectCandidateDplus = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagHadron;
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagHadron || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagHadron);

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry trackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext&)
  {
    std::array<bool, 20> processes = {doprocessDataDplusToPiKPi, doprocessMcDplusToPiKPi, doprocessDataDplusToPiKPiWithML, doprocessMcDplusToPiKPiWithML, doprocessMcDplusToPiKPiGen,
                                      doprocessDataLcToPKPi, doprocessMcLcToPKPi, doprocessDataLcToPKPiWithML, doprocessMcLcToPKPiWithML, doprocessMcLcToPKPiGen, doprocessDataD0ToPiK, doprocessMcD0ToPiK, doprocessDataD0ToPiKWithML, doprocessMcD0ToPiKWithML, doprocessMcD0ToPiKGen, doprocessDataDstarToD0Pi, doprocessMcDstarToD0Pi, doprocessDataDstarToD0PiWithML, doprocessMcDstarToD0PiWithML, doprocessMcDstarToD0PiGen};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int const cutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    trackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});

    // event QA histograms
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
    magField = 0.0;
    /// Initializing CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    hfEvSel.init(qaRegistry, &zorroSummary); // collision monitoring

    int64_t const now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    bool useLcMl = doprocessDataLcToPKPiWithML || doprocessMcLcToPKPiWithML;
    bool useDplusMl = doprocessDataDplusToPiKPiWithML || doprocessMcDplusToPiKPiWithML;
    bool useD0Ml = doprocessDataD0ToPiKWithML || doprocessMcD0ToPiKWithML;
    bool useDstarMl = doprocessDataDstarToD0PiWithML || doprocessMcDstarToD0PiWithML;

    if (applyMlMode == FillMlFromNewBDT) {
      if (useLcMl) {
        hfMlResponseLc.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        hfMlResponseLc.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseLc.init();
      }
      if (useDplusMl) {
        hfMlResponseDplus.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        hfMlResponseDplus.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseDplus.init();
      }
      if (useD0Ml) {
        hfMlResponseD0.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        hfMlResponseD0.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseD0.init();
      }
      if (useDstarMl) {
        hfMlResponseDstar.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        hfMlResponseDstar.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseDstar.init();
      }

      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        if (useLcMl) {
          hfMlResponseLc.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        }
        if (useDplusMl) {
          hfMlResponseDplus.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        }
        if (useD0Ml) {
          hfMlResponseD0.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        }
        if (useDstarMl) {
          hfMlResponseDstar.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        }

      } else {
        if (useLcMl) {
          hfMlResponseLc.setModelPathsLocal(onnxFileNames);
        }
        if (useDplusMl) {
          hfMlResponseDplus.setModelPathsLocal(onnxFileNames);
        }
        if (useD0Ml) {
          hfMlResponseD0.setModelPathsLocal(onnxFileNames);
        }
        if (useDstarMl) {
          hfMlResponseDstar.setModelPathsLocal(onnxFileNames);
        }
      }
    }
  }

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(const aod::BCsWithTimestamps::iterator& bc)
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
                     -999., -999., -999., -999.,
                     -999., -999., -999., -999.,
                     -999., -999., -999.);
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
      constexpr int GenFromTransport = -1; // -1 if a particle produced during transport
      // get list of mothers, but it could be empty (for example in case of injected light nuclei)
      auto motherparticlesMc = particleMc.template mothers_as<aod::McParticles>();
      // check pdg code
      // if this fails, the particle is a fake
      if (std::abs(pdgCode) == std::abs(trkPDGCode.value)) {
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
        } else if (particleMc.getProcess() == TMCProcess::kPDecay && particleMc.getGenStatusCode() == GenFromTransport && !motherparticlesMc.empty()) {
          // get direct mother
          auto motherparticleMc = motherparticlesMc.front();
          pdgCodeMother = motherparticleMc.pdgCode();
          particleOrigin = checkDaughterType(fdparttype, motherparticleMc.pdgCode());
          // check if particle is material
          // particle is from inelastic hadronic interaction -> getProcess() == 23
          // particle is generated during transport -> getGenStatusCode() == -1
        } else if (particleMc.getProcess() == TMCProcess::kPHInhelastic && particleMc.getGenStatusCode() == GenFromTransport) {
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

  template <bool IsMc = false, typename TrackType, typename CollisionType>
  bool fillTracksForCharmHadron(CollisionType const& col, TrackType const& tracks)
  {

    std::vector<int> const childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    // std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    bool fIsTrackFilled = false;

    for (const auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, true>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<false, aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));
      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      int64_t timeStamp = bc.timestamp();
      // track global index
      outputPartsIndex(track.globalIndex());
      outputPartsTime(timeStamp);
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

      if constexpr (IsMc) {
        fillMcParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }
    return fIsTrackFilled;
  }

  template <DecayChannel Channel, bool IsMc, bool UseCharmMl, typename TrackType, typename CollisionType, typename CandType>
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

    qaRegistry.fill(HIST("hEventQA"), 1 + Event::All);

    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(col, rejectionMask, mult);
    if (rejectionMask != 0) {
      /// at least one event selection not satisfied --> reject the candidate
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

    // Filling candidate properties
    if constexpr (Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::LcToPKPi) {
      rowCandCharm3Prong.reserve(sizeCand);
    } else if constexpr (Channel == DecayChannel::D0ToPiK) {
      rowCandCharm2Prong.reserve(sizeCand);
    } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
      rowCandCharmDstar.reserve(sizeCand);
    }
    bool isTrackFilled = false;
    bool isSelectedMlLcToPKPi = true;
    bool isSelectedMlLcToPiKP = true;
    bool isSelectedMlDplusToPiKPi = true;
    bool isSelectedMlD0ToPiK = true;
    bool isSelectedMlD0barToKPi = true;
    bool isSelectedMlDstarToD0Pi = true;

    for (const auto& candidate : candidates) {
      outputMlD0 = {-1.0f, -1.0f, -1.0f};
      outputMlD0bar = {-1.0f, -1.0f, -1.0f};
      outputMlDplus = {-1.0f, -1.0f, -1.0f};
      outputMlDstar = {-1.0f, -1.0f, -1.0f};
      outputMlPKPi = {-1.0f, -1.0f, -1.0f};
      outputMlPiKP = {-1.0f, -1.0f, -1.0f};
      auto trackPos1 = candidate.template prong0_as<TrackType>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<TrackType>();  // negative daughter (positive for the antiparticles)

      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      int64_t timeStamp = bc.timestamp();

      auto fillTable = [&](int candFlag,
                           int functionSelection,
                           float bdtScoreBkg,
                           float bdtScorePrompt,
                           float bdtScoreFd) {
        if (functionSelection >= 1) {
          if constexpr (Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::LcToPKPi) {
            auto trackPos2 = candidate.template prong2_as<TrackType>();
            rowCandCharm3Prong(
              outputCollision.lastIndex(),
              timeStamp,
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
              1 << candFlag,
              bdtScoreBkg,
              bdtScorePrompt,
              bdtScoreFd);

          } else if constexpr (Channel == DecayChannel::D0ToPiK) {
            int signD0 = -999;
            if (candFlag == static_cast<int>(D0CandFlag::D0)) {
              signD0 = +1;
            } else if (candFlag == static_cast<int>(D0CandFlag::D0Bar)) {
              signD0 = -1;
            } else if (candFlag == static_cast<int>(D0CandFlag::Reflected)) {
              signD0 = 0;
            } else {
              LOG(error) << "Unexpected candFlag = " << candFlag;
            }
            rowCandCharm2Prong(
              outputCollision.lastIndex(),
              timeStamp,
              signD0,
              trackPos1.globalIndex(),
              trackNeg.globalIndex(),
              trackPos1.pt(),
              trackNeg.pt(),
              trackPos1.eta(),
              trackNeg.eta(),
              trackPos1.phi(),
              trackNeg.phi(),
              1 << candFlag,
              bdtScoreBkg,
              bdtScorePrompt,
              bdtScoreFd);
          } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
            auto trackPos2 = candidate.template prongPi_as<TrackType>();
            rowCandCharmDstar(
              outputCollision.lastIndex(),
              timeStamp,
              candidate.template prongPi_as<TrackType>().sign(),
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
              1 << candFlag,
              bdtScoreBkg,
              bdtScorePrompt,
              bdtScoreFd);
          }

          if constexpr (IsMc) {
            rowCandMcCharmHad(
              candidate.flagMcMatchRec(),
              candidate.originMcRec());
          }
        }
      };

      if constexpr (Channel == DecayChannel::DplusToPiKPi) {
        if constexpr (UseCharmMl) {
          /// fill with ML information
          /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
          if (applyMlMode == FillMlFromSelector) {
            if (candidate.mlProbDplusToPiKPi().size() > 0) {
              outputMlDplus.at(0) = candidate.mlProbDplusToPiKPi()[0]; /// bkg score
              outputMlDplus.at(1) = candidate.mlProbDplusToPiKPi()[1]; /// prompt score
              outputMlDplus.at(2) = candidate.mlProbDplusToPiKPi()[2]; /// non-prompt score
            }
          } else if (applyMlMode == FillMlFromNewBDT) {
            isSelectedMlDplusToPiKPi = false;
            if (candidate.mlProbDplusToPiKPi().size() > 0) {
              std::vector<float> inputFeaturesDplusToPiKPi = hfMlResponseDplus.getInputFeatures(candidate);
              isSelectedMlDplusToPiKPi = hfMlResponseDplus.isSelectedMl(inputFeaturesDplusToPiKPi, candidate.pt(), outputMlDplus);
            }
            if (!isSelectedMlDplusToPiKPi) {
              continue;
            }
          } else {
            LOGF(fatal, "Please check your Ml configuration!!");
          }
        }
        fillTable(2, candidate.isSelDplusToPiKPi(), outputMlDplus.at(0), outputMlDplus.at(1), outputMlDplus.at(2));

      } else if constexpr (Channel == DecayChannel::LcToPKPi) {
        if constexpr (UseCharmMl) {
          /// fill with ML information
          /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
          if (applyMlMode == FillMlFromSelector) {
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
          } else if (applyMlMode == FillMlFromNewBDT) {
            isSelectedMlLcToPKPi = false;
            isSelectedMlLcToPiKP = false;
            if (candidate.mlProbLcToPKPi().size() > 0) {
              std::vector<float> inputFeaturesLcToPKPi = hfMlResponseLc.getInputFeatures(candidate, true);
              isSelectedMlLcToPKPi = hfMlResponseLc.isSelectedMl(inputFeaturesLcToPKPi, candidate.pt(), outputMlPKPi);
            }
            if (candidate.mlProbLcToPiKP().size() > 0) {
              std::vector<float> inputFeaturesLcToPiKP = hfMlResponseLc.getInputFeatures(candidate, false);
              isSelectedMlLcToPiKP = hfMlResponseLc.isSelectedMl(inputFeaturesLcToPiKP, candidate.pt(), outputMlPKPi);
            }
            if (!isSelectedMlLcToPKPi && !isSelectedMlLcToPiKP) {
              continue;
            }
          } else {
            LOGF(fatal, "Please check your Ml configuration!!");
          }
        }
        fillTable(0, candidate.isSelLcToPKPi(), outputMlPKPi.at(0), outputMlPKPi.at(1), outputMlPKPi.at(2));
        fillTable(1, candidate.isSelLcToPiKP(), outputMlPiKP.at(0), outputMlPiKP.at(1), outputMlPiKP.at(2));
      } else if constexpr (Channel == DecayChannel::D0ToPiK) {
        if constexpr (UseCharmMl) {

          /// fill with ML information
          /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
          if (applyMlMode == FillMlFromSelector) {
            if (candidate.mlProbD0().size() > 0) {
              outputMlD0.at(0) = candidate.mlProbD0()[0]; /// bkg score
              outputMlD0.at(1) = candidate.mlProbD0()[1]; /// prompt score
              outputMlD0.at(2) = candidate.mlProbD0()[2]; /// non-prompt score
            }
            if (candidate.mlProbD0bar().size() > 0) {
              outputMlD0bar.at(0) = candidate.mlProbD0bar()[0]; /// bkg score
              outputMlD0bar.at(1) = candidate.mlProbD0bar()[1]; /// prompt score
              outputMlD0bar.at(2) = candidate.mlProbD0bar()[2]; /// non-prompt score
            }

          } else if (applyMlMode == FillMlFromNewBDT) {
            isSelectedMlD0ToPiK = false;
            isSelectedMlD0barToKPi = false;

            if (candidate.mlProbD0().size() > 0) {
              std::vector<float> inputFeaturesD0ToPiK = hfMlResponseD0.getInputFeatures(candidate, o2::constants::physics::kD0);
              isSelectedMlD0ToPiK = hfMlResponseD0.isSelectedMl(inputFeaturesD0ToPiK, candidate.pt(), outputMlD0);
            }
            if (candidate.mlProbD0bar().size() > 0) {
              std::vector<float> inputFeaturesD0barToKPi = hfMlResponseD0.getInputFeatures(candidate, o2::constants::physics::kD0);
              isSelectedMlD0barToKPi = hfMlResponseD0.isSelectedMl(inputFeaturesD0barToKPi, candidate.pt(), outputMlD0bar);
            }
            if (!isSelectedMlD0ToPiK && !isSelectedMlD0barToKPi) {
              continue;
            }
          } else {
            LOGF(fatal, "Please check your Ml configuration!!");
          }
        }
        if (candidate.isSelD0() && candidate.isSelD0bar()) {
          fillTable(2, candidate.isSelD0(), outputMlD0.at(0), outputMlD0.at(1), outputMlD0.at(2)); // tag reflection
        } else {
          fillTable(0, candidate.isSelD0(), outputMlD0.at(0), outputMlD0.at(1), outputMlD0.at(2));
          fillTable(1, candidate.isSelD0bar(), outputMlD0bar.at(0), outputMlD0bar.at(1), outputMlD0bar.at(2));
        }

      } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
        if constexpr (UseCharmMl) {
          /// fill with ML information
          /// BDT index 0: bkg score; BDT index 1: prompt score; BDT index 2: non-prompt score
          if (applyMlMode == FillMlFromSelector) {
            if (candidate.mlProbDstarToD0Pi().size() > 0) {
              outputMlDstar.at(0) = candidate.mlProbDstarToD0Pi()[0]; /// bkg score
              outputMlDstar.at(1) = candidate.mlProbDstarToD0Pi()[1]; /// prompt score
              outputMlDstar.at(2) = candidate.mlProbDstarToD0Pi()[2]; /// non-prompt score
            }
          } else if (applyMlMode == FillMlFromNewBDT) {
            isSelectedMlDstarToD0Pi = false;
            if (candidate.mlProbDstarToD0Pi().size() > 0) {
              std::vector<float> inputFeaturesDstarToD0Pi = hfMlResponseDstar.getInputFeatures(candidate, false);
              isSelectedMlDstarToD0Pi = hfMlResponseDstar.isSelectedMl(inputFeaturesDstarToD0Pi, candidate.pt(), outputMlDstar);
            }
            if (!isSelectedMlDstarToD0Pi) {
              continue;
            }
          } else {
            LOGF(fatal, "Please check your Ml configuration!!");
          }
        }
        fillTable(2, candidate.isSelDstarToD0Pi(), outputMlDstar.at(0), outputMlDstar.at(1), outputMlDstar.at(2));
      }
    }
    isTrackFilled = fillTracksForCharmHadron<IsMc>(col, tracks);

    aod::femtodreamcollision::BitMaskType bitTrack = 0;
    if (isTrackFilled) {
      bitTrack |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::TrackSelected);
    }

    aod::femtodreamcollision::BitMaskType bitCand = 0;
    if (sizeCand > 0) {
      bitCand |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::CharmSelected);
    }

    if (isTrackFilled && (sizeCand > 0)) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::PairSelected);
    }

    rowMasks(bitTrack,
             bitCand,
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

  template <DecayChannel Channel, typename ParticleType>
  void fillCharmHadMcGen(ParticleType particles)
  {
    // Filling particle properties
    rowCandCharmHadGen.reserve(particles.size());
    if constexpr (Channel == DecayChannel::DplusToPiKPi) {
      for (const auto& particle : particles) {
        if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
          rowCandCharmHadGen(
            particle.mcCollisionId(),
            particle.flagMcMatchGen(),
            particle.originMcGen());
        }
      }
    } else if constexpr (Channel == DecayChannel::LcToPKPi) {
      for (const auto& particle : particles) {
        if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
          rowCandCharmHadGen(
            particle.mcCollisionId(),
            particle.flagMcMatchGen(),
            particle.originMcGen());
        }
      }
    } else if constexpr (Channel == DecayChannel::D0ToPiK) {
      for (const auto& particle : particles) {
        if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          rowCandCharmHadGen(
            particle.mcCollisionId(),
            particle.flagMcMatchGen(),
            particle.originMcGen());
        }
      }
    } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
      for (const auto& particle : particles) {
        if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) {
          rowCandCharmHadGen(
            particle.mcCollisionId(),
            particle.flagMcMatchGen(),
            particle.originMcGen());
        }
      }
    }
  }

  /// D0ToPiK
  void processDataD0ToPiK(FemtoFullCollision const& col,
                          aod::BCsWithTimestamps const&,
                          FemtoHFTracks const& tracks,
                          soa::Filtered<CandidateD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataD0ToPiK, "Provide experimental data for D0ToPiK femto", false);

  void processDataD0ToPiKWithML(FemtoFullCollision const& col,
                                aod::BCsWithTimestamps const&,
                                FemtoHFTracks const& tracks,
                                soa::Filtered<soa::Join<CandidateD0,
                                                        aod::HfMlD0>> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataD0ToPiKWithML, "Provide experimental data for D0ToPiK with ml", false);

  void processMcD0ToPiK(FemtoFullCollisionMc const& col,
                        aod::BCsWithTimestamps const&,
                        FemtoHFMcTracks const& tracks,
                        aod::McParticles const&,
                        CandidateD0Mc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcD0ToPiK, "Provide Mc for D0ToPiK", false);

  void processMcD0ToPiKWithML(FemtoFullCollisionMc const& col,
                              aod::BCsWithTimestamps const&,
                              FemtoHFMcTracks const& tracks,
                              aod::McParticles const&,
                              soa::Join<CandidateD0Mc,
                                        aod::HfMlD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcD0ToPiKWithML, "Provide Mc for D0ToPiK with ml", false);

  void processMcD0ToPiKGen(Generated2ProngMc const& particles)
  {
    fillCharmHadMcGen<DecayChannel::D0ToPiK>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcD0ToPiKGen, "Provide Mc Generated D0ToPiK", false);

  /// DstarToD0Pi
  void processDataDstarToD0Pi(FemtoFullCollision const& col,
                              aod::BCsWithTimestamps const&,
                              FemtoHFTracks const& tracks,
                              soa::Filtered<CandidateDstar> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataDstarToD0Pi, "Provide experimental data for DstarToD0Pi femto", false);

  void processDataDstarToD0PiWithML(FemtoFullCollision const& col,
                                    aod::BCsWithTimestamps const&,
                                    FemtoHFTracks const& tracks,
                                    soa::Filtered<soa::Join<CandidateDstar,
                                                            aod::HfMlDstarToD0Pi>> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataDstarToD0PiWithML, "Provide experimental data for DstarToD0Pi with ml", false);

  void processMcDstarToD0Pi(FemtoFullCollisionMc const& col,
                            aod::BCsWithTimestamps const&,
                            FemtoHFMcTracks const& tracks,
                            aod::McParticles const&,
                            CandidateDstarMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDstarToD0Pi, "Provide Mc for DstarToD0Pi", false);

  void processMcDstarToD0PiWithML(FemtoFullCollisionMc const& col,
                                  aod::BCsWithTimestamps const&,
                                  FemtoHFMcTracks const& tracks,
                                  aod::McParticles const&,
                                  soa::Join<CandidateDstarMc,
                                            aod::HfMlDstarToD0Pi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDstarToD0PiWithML, "Provide Mc for DstarToD0Pi with ml", false);

  void processMcDstarToD0PiGen(GeneratedDstarMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::DstarToD0Pi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDstarToD0PiGen, "Provide Mc Generated DstarToD0Pi", false);

  /// DplusToPiKPi
  void processDataDplusToPiKPi(FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               FemtoHFTracks const& tracks,
                               soa::Filtered<CandidateDplus> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataDplusToPiKPi, "Provide experimental data for DplusToPiKPi femto", false);

  void processDataDplusToPiKPiWithML(FemtoFullCollision const& col,
                                     aod::BCsWithTimestamps const&,
                                     FemtoHFTracks const& tracks,
                                     soa::Filtered<soa::Join<CandidateDplus,
                                                             aod::HfMlDplusToPiKPi>> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataDplusToPiKPiWithML, "Provide experimental data for DplusToPiKPi with ml", false);

  void processMcDplusToPiKPi(FemtoFullCollisionMc const& col,
                             aod::BCsWithTimestamps const&,
                             FemtoHFMcTracks const& tracks,
                             aod::McParticles const&,
                             CandidateDplusMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDplusToPiKPi, "Provide Mc for DplusToPiKPi", false);

  void processMcDplusToPiKPiWithML(FemtoFullCollisionMc const& col,
                                   aod::BCsWithTimestamps const&,
                                   FemtoHFMcTracks const& tracks,
                                   aod::McParticles const&,
                                   soa::Join<CandidateDplusMc,
                                             aod::HfMlDplusToPiKPi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDplusToPiKPiWithML, "Provide Mc for DplusToPiKPi with ml", false);

  void processMcDplusToPiKPiGen(Generated3ProngMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::DplusToPiKPi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcDplusToPiKPiGen, "Provide Mc Generated DplusToPiKPi", false);

  /// LcToPKPi
  void processDataLcToPKPi(FemtoFullCollision const& col,
                           aod::BCsWithTimestamps const&,
                           FemtoHFTracks const& tracks,
                           soa::Filtered<CandidateLc> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, false, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataLcToPKPi, "Provide experimental data for Lc(PKPi)-proton femto", false);

  void processDataLcToPKPiWithML(FemtoFullCollision const& col,
                                 aod::BCsWithTimestamps const&,
                                 FemtoHFTracks const& tracks,
                                 soa::Filtered<soa::Join<CandidateLc,
                                                         aod::HfMlLcToPKPi>> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, false, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processDataLcToPKPiWithML, "Provide experimental data for Lc(PKPi)-proton femto with ml", false);

  void processMcLcToPKPi(FemtoFullCollisionMc const& col,
                         aod::BCsWithTimestamps const&,
                         FemtoHFMcTracks const& tracks,
                         aod::McParticles const&,
                         CandidateLcMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, true, false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcLcToPKPi, "Provide Mc for lctopkpi", false);

  void processMcLcToPKPiWithML(FemtoFullCollisionMc const& col,
                               aod::BCsWithTimestamps const&,
                               FemtoHFMcTracks const& tracks,
                               aod::McParticles const&,
                               soa::Join<CandidateLcMc,
                                         aod::HfMlLcToPKPi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, true, true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcLcToPKPiWithML, "Provide Mc for lctopkpi with ml", false);

  void processMcLcToPKPiGen(Generated3ProngMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::LcToPKPi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsTrackFemtoDream, processMcLcToPKPiGen, "Provide Mc Generated lctopkpi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfProducerCharmHadronsTrackFemtoDream>(cfgc)};
}
