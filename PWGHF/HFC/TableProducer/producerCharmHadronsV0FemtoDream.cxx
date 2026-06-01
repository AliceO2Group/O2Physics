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

/// \file producerCharmHadronsV0FemtoDream.cxx
/// \brief Tasks that produces the V0 and CharmHadrons tables used for the pairing
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGCF/FemtoDream/Core/femtoDreamV0SelectionK0Short.h"
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
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
#include <cstddef>
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
  RejNoV0sAndCharm,
  V0Selected,
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

enum V0Channel {
  None = 0,
  K0S,
  Lambda
};

enum class D0CandFlag : uint8_t {
  D0 = 0,
  D0Bar = 1,
  Reflected = 2
};

struct HfProducerCharmHadronsV0FemtoDream {

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
  Configurable<int> selectionFlagV0{"selectionFlagV0", 1, "Activate filling of V0: 0 for none, 1 for K0S,  2 for Lambda"};
  /// Charm hadron table
  Configurable<int> selectionFlagCharmHadron{"selectionFlagCharmHadron", 1, "Selection Flag for Charm Hadron: 1 for Lc, 7 for Dplus (Topologic and PID cuts)"};
  Configurable<bool> useCent{"useCent", false, "Enable centrality for Charm Hadron"};

  Configurable<int> v0PDGCode{"v0PDGCode", 310, "PDG code of the selected V0 (310: K0S, 3122: Lambda) for Monte Carlo truth "};
  Configurable<float> trkPIDnSigmaOffsetTPC{"trkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"}; // set to zero for run3 or so
  Configurable<float> trkPIDnSigmaOffsetTOF{"trkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<float> trkPIDnSigmaForLambdaSign{"trkPIDnSigmaForLambdaSign", 0., "daught track PID for Lambda sign determination"};

  Configurable<bool> trkRejectNotPropagated{"trkRejectNotPropagated", false, "True: reject not propagated tracks"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confLambdaSign{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0Sign, "confLambda"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaPtMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMin, "confLambda"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaPtMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMax, "confLambda"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaEtaMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0etaMax, "confLambda"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaDCADaughMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DCADaughMax, "confLambda"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaCPAMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0CPAMin, "confLambda"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaTranRadMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMin, "confLambda"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaTranRadMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMax, "confLambda"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> confLambdaDecVtxMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DecVtxMax, "confLambda"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<float> confLambdaInvMassLowLimit{"confLambdaInvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> confLambdaInvMassUpLimit{"confLambdaInvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};
    Configurable<bool> confLambdaRejectKaons{"confLambdaRejectKaons", false, "Switch to reject kaons"};
    Configurable<bool> confLambdaRejectLambdas{"confLambdaRejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
    Configurable<float> confLambdaInvKaonMassLowLimit{"confLambdaInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> confLambdaInvKaonMassUpLimit{"confLambdaInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    Configurable<std::vector<float>> confLambdaChildSign{"confLambdaChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> confLambdaChildEtaMax{"confLambdaChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> confLambdaChildTPCnClsMin{"confLambdaChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confLambdaChildDCAMin{"confLambdaChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confLambdaChildPIDnSigmaMax{"confLambdaChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confLambdaChildPIDspecies{"confLambdaChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

    // cuts and object for v0 2
    Configurable<std::vector<float>> confK0shortSign{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0Sign, "confK0short"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortPtMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMin, "confK0short"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortPtMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMax, "confK0short"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortEtaMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0etaMax, "confK0short"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortDCADaughMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DCADaughMax, "confK0short"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortCPAMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0CPAMin, "confK0short"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortTranRadMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMin, "confK0short"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortTranRadMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMax, "confK0short"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> confK0shortDecVtxMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DecVtxMax, "confK0short"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<float> confK0shortInvMassLowLimit{"confK0shortInvMassLowLimit", 0.45, "Lower limit of the V0 invariant mass"};
    Configurable<float> confK0shortInvMassUpLimit{"confK0shortInvMassUpLimit", 0.55, "Upper limit of the V0 invariant mass"};
    Configurable<bool> confK0shortRejectKaons{"confK0shortRejectKaons", false, "Switch to reject kaons"};
    Configurable<bool> confK0shortRejectLambdas{"confK0shortRejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
    Configurable<float> confK0shortInvKaonMassLowLimit{"confK0shortInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> confK0shortInvKaonMassUpLimit{"confK0shortInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    Configurable<std::vector<float>> confK0shortChildSign{"confK0shortChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> confK0shortChildEtaMax{"confK0shortChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> confK0shortChildTPCnClsMin{"confK0shortChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confK0shortChildDCAMin{"confK0shortChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confK0shortChildPIDnSigmaMax{"confK0shortChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confK0shortChildPIDspecies{"confK0shortChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Pion}, "V0 Child sel: Particles species for PID"};
  } V0Sel;

  // ML inference
  Configurable<int> applyMlMode{"applyMlMode", 1, "None: 0, BDT model from candidate selector: 1, New BDT model on Top of candidate selector: 2"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

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
  using FemtoHFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullEl, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullEl, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using FemtoHFMcTracks = soa::Join<aod::McTrackLabels, FemtoHFTracks>;
  using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
  using FemtoFullCollisionMc = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
  using FemtoFullMcgenCollisions = soa::Join<aod::McCollisions, o2::aod::MultsExtraMC>;
  using FemtoFullMcgenCollision = FemtoFullMcgenCollisions::iterator;

  using Generated3ProngMc = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using Generated2ProngMc = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  using GeneratedDstarMc = soa::Join<aod::McParticles, aod::HfCandDstarMcGen>;

  Filter filterSelectCandidateD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagCharmHadron || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagCharmHadron);
  Filter filterSelectCandidateDstar = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;
  Filter filterSelectCandidateDplus = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagCharmHadron;
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagCharmHadron || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagCharmHadron);

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistryV0{"V0", {}, OutputObjHandlingPolicy::AnalysisObject};
  FemtoDreamV0Selection k0sCuts;
  FemtoDreamV0Selection lambdaCuts;
  void init(InitContext&)
  {
    std::array<bool, 20> processes = {doprocessDataDplusToPiKPi, doprocessMcDplusToPiKPi, doprocessDataDplusToPiKPiWithML, doprocessMcDplusToPiKPiWithML, doprocessMcDplusToPiKPiGen,
                                      doprocessDataLcToPKPi, doprocessMcLcToPKPi, doprocessDataLcToPKPiWithML, doprocessMcLcToPKPiWithML, doprocessMcLcToPKPiGen, doprocessDataD0ToPiK, doprocessMcD0ToPiK, doprocessDataD0ToPiKWithML, doprocessMcD0ToPiKWithML, doprocessMcD0ToPiKGen, doprocessDataDstarToD0Pi, doprocessMcDstarToD0Pi, doprocessDataDstarToD0PiWithML, doprocessMcDstarToD0PiWithML, doprocessMcDstarToD0PiGen};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int const cutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    qaRegistryV0.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});

    // event QA histograms
    constexpr int EventTypes = PairSelected + 1;
    std::string labels[EventTypes];
    labels[Event::All] = "All events";
    labels[Event::RejEveSel] = "rejected by event selection";
    labels[Event::RejNoV0sAndCharm] = "rejected by no V0s and charm";
    labels[Event::V0Selected] = "with V0 ";
    labels[Event::CharmSelected] = "with charm hadrons ";
    labels[Event::PairSelected] = "with pairs";

    static const AxisSpec axisEvents = {EventTypes, 0.5, EventTypes + 0.5, ""};
    qaRegistry.add("hEventQA", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < EventTypes; iBin++) {
      qaRegistry.get<TH1>(HIST("hEventQA"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    if (selectionFlagV0 == V0Channel::K0S) {
      k0sCuts.setSelection(V0Sel.confK0shortSign, femto_dream_v0_selection::kV0Sign, femtoDreamSelection::kEqual);
      k0sCuts.setSelection(V0Sel.confK0shortPtMin, femto_dream_v0_selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      k0sCuts.setSelection(V0Sel.confK0shortPtMax, femto_dream_v0_selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      k0sCuts.setSelection(V0Sel.confK0shortEtaMax, femto_dream_v0_selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      k0sCuts.setSelection(V0Sel.confK0shortDCADaughMax, femto_dream_v0_selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      k0sCuts.setSelection(V0Sel.confK0shortCPAMin, femto_dream_v0_selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      k0sCuts.setSelection(V0Sel.confK0shortTranRadMin, femto_dream_v0_selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      k0sCuts.setSelection(V0Sel.confK0shortTranRadMax, femto_dream_v0_selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      k0sCuts.setSelection(V0Sel.confK0shortDecVtxMax, femto_dream_v0_selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      k0sCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      k0sCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      k0sCuts.setChildPIDSpecies(femto_dream_v0_selection::kPosTrack, V0Sel.confK0shortChildPIDspecies);
      k0sCuts.setChildPIDSpecies(femto_dream_v0_selection::kNegTrack, V0Sel.confK0shortChildPIDspecies);
      k0sCuts.init<aod::femtodreamparticle::ParticleType::kV0K0Short, aod::femtodreamparticle::ParticleType::kV0K0ShortChild, aod::femtodreamparticle::cutContainerType>(&qaRegistryV0, &qaRegistryV0);
      k0sCuts.setInvMassLimits(V0Sel.confK0shortInvMassLowLimit, V0Sel.confK0shortInvMassUpLimit);
      k0sCuts.setIsMother(false);

      k0sCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kPosTrack, trkRejectNotPropagated);
      k0sCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kNegTrack, trkRejectNotPropagated);

      k0sCuts.setnSigmaPIDOffsetTPC(trkPIDnSigmaOffsetTPC);
      k0sCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kPosTrack, trkPIDnSigmaOffsetTPC, trkPIDnSigmaOffsetTOF);
      k0sCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kNegTrack, trkPIDnSigmaOffsetTPC, trkPIDnSigmaOffsetTOF);

      k0sCuts.setRejectLambda(V0Sel.confK0shortRejectLambdas);
      k0sCuts.setKaonInvMassLimits(V0Sel.confK0shortInvKaonMassLowLimit, V0Sel.confK0shortInvKaonMassUpLimit);
    }

    if (selectionFlagV0 == V0Channel::Lambda) {
      // lambdaCuts.setSelection(V0Sel.confLambdaSign, femto_dream_v0_selection::kV0Sign, femtoDreamSelection::kEqual);
      lambdaCuts.setSelection(V0Sel.confLambdaPtMin, femto_dream_v0_selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaPtMax, femto_dream_v0_selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaEtaMax, femto_dream_v0_selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaDCADaughMax, femto_dream_v0_selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaCPAMin, femto_dream_v0_selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaTranRadMin, femto_dream_v0_selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaTranRadMax, femto_dream_v0_selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      lambdaCuts.setSelection(V0Sel.confLambdaDecVtxMax, femto_dream_v0_selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      lambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      lambdaCuts.setChildCuts(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      lambdaCuts.setChildPIDSpecies(femto_dream_v0_selection::kPosTrack, V0Sel.confLambdaChildPIDspecies);
      lambdaCuts.setChildPIDSpecies(femto_dream_v0_selection::kNegTrack, V0Sel.confLambdaChildPIDspecies);
      lambdaCuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistryV0, &qaRegistryV0);
      lambdaCuts.setInvMassLimits(V0Sel.confLambdaInvMassLowLimit, V0Sel.confLambdaInvMassUpLimit);
      lambdaCuts.setIsMother(true);

      lambdaCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kPosTrack, trkRejectNotPropagated);
      lambdaCuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kNegTrack, trkRejectNotPropagated);

      lambdaCuts.setnSigmaPIDOffsetTPC(trkPIDnSigmaOffsetTPC);
      lambdaCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kPosTrack, trkPIDnSigmaOffsetTPC, trkPIDnSigmaOffsetTOF);
      lambdaCuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kNegTrack, trkPIDnSigmaOffsetTPC, trkPIDnSigmaOffsetTOF);

      if (V0Sel.confLambdaRejectKaons) {
        lambdaCuts.setKaonInvMassLimits(V0Sel.confLambdaInvKaonMassLowLimit, V0Sel.confLambdaInvKaonMassUpLimit);
      }
    }

    runNumber = 0;
    magField = 0.0;
    /// Initializing CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    hfEvSel.addHistograms(qaRegistry); // collision monitoring

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

  template <bool isTrackOrV0, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle, int const& signV0 = 0)
  {
    if constexpr (isTrackOrV0) {
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
    } else {
      auto sign = signV0;
      outputDebugParts(sign,                                                          // sign
                       -999., -999., -999., -999., -999., -999., -999., -999., -999., // track properties (DCA, NCls, crossed rows, etc.)
                       -999., -999., -999., -999., -999., -999., -999., -999.,        // TPC PID (TPC signal + particle hypothesis)
                       -999., -999., -999., -999., -999., -999., -999.,               // TOF PID
                       -999., -999., -999., -999., -999., -999., -999., -999.,        // ITS PID
                       particle.dcaV0daughters(),
                       particle.v0radius(),
                       particle.x(),
                       particle.y(),
                       particle.z(),
                       particle.mK0Short(),
                       -999., -999., -999., -999., -999., -999., -999.);
    }
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
      if (std::abs(pdgCode) == std::abs(v0PDGCode.value)) {
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

  template <bool IsMc = false, typename V0Type, typename CollisionType, typename TrackType>
  bool fillV0sForCharmHadron(CollisionType const& col, V0Type const& fullV0s, TrackType const& /*tracks*/)
  {
    const bool isK0S = (selectionFlagV0 == V0Channel::K0S);
    const bool isLambda = (selectionFlagV0 == V0Channel::Lambda);
    if (!isK0S && !isLambda) {
      LOG(fatal) << "Invalid V0 particle !! Please check the configuration";
    }

    std::vector<int> childIDs = {0, 0};
    std::vector<int> tmpIDtrack;
    bool fIsV0Filled = false;

    auto bc = col.template bc_as<aod::BCsWithTimestamps>();
    int64_t timeStamp = bc.timestamp();

    for (const auto& v0 : fullV0s) {

      auto postrack = v0.template posTrack_as<TrackType>();
      auto negtrack = v0.template negTrack_as<TrackType>();

      float massV0 = 0.f;
      float antiMassV0 = 0.f;

      int signV0 = determineV0Sign(v0, postrack, negtrack);

      std::array<aod::femtodreamparticle::cutContainerType, 5> cutContainerV0{};

      if (isK0S) {
        k0sCuts.fillLambdaQA<aod::femtodreamparticle::ParticleType::kV0K0Short>(col, v0, postrack, negtrack);
        if (!k0sCuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }
        k0sCuts.fillQA<aod::femtodreamparticle::ParticleType::kV0K0Short,
                       aod::femtodreamparticle::ParticleType::kV0K0ShortChild>(col, v0, postrack, negtrack);

        cutContainerV0 = k0sCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);
        massV0 = v0.mK0Short();
        antiMassV0 = v0.mK0Short();
      } else { // Lambda
        lambdaCuts.fillLambdaQA<aod::femtodreamparticle::ParticleType::kV0>(col, v0, postrack, negtrack);
        if (!lambdaCuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }
        lambdaCuts.fillQA<aod::femtodreamparticle::ParticleType::kV0,
                          aod::femtodreamparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack);
        cutContainerV0 = lambdaCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);
        massV0 = v0.mLambda();
        antiMassV0 = v0.mAntiLambda();
      }

      bool passPosPID = false, passNegPID = false;
      isV0DaughterPidSelected(postrack, negtrack, signV0, passPosPID, passNegPID);

      // --- pos child
      int rowPos = getRowDaughters(v0.posTrackId(), tmpIDtrack);
      childIDs[0] = rowPos;
      childIDs[1] = 0;

      auto daughType = isK0S ? aod::femtodreamparticle::ParticleType::kV0K0ShortChild
                             : aod::femtodreamparticle::ParticleType::kV0Child;
      outputPartsTime(timeStamp);
      outputPartsIndex(v0.posTrackId());
      outputParts(outputCollision.lastIndex(),
                  v0.positivept(), v0.positiveeta(), v0.positivephi(),
                  daughType,
                  cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kPosCuts),
                  static_cast<aod::femtodreamparticle::cutContainerType>(passPosPID),
                  postrack.dcaXY(),
                  childIDs,
                  0,
                  0);
      const int rowOfPosTrack = outputParts.lastIndex();

      // --- neg child
      int rowNeg = getRowDaughters(v0.negTrackId(), tmpIDtrack);
      childIDs[0] = 0;
      childIDs[1] = rowNeg;
      outputPartsTime(timeStamp);
      outputPartsIndex(v0.negTrackId());
      outputParts(outputCollision.lastIndex(),
                  v0.negativept(), v0.negativeeta(), v0.negativephi(),
                  daughType,
                  cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kNegCuts),
                  static_cast<aod::femtodreamparticle::cutContainerType>(passNegPID),
                  negtrack.dcaXY(),
                  childIDs,
                  0,
                  0);
      const int rowOfNegTrack = outputParts.lastIndex();
      // --- mother
      std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};
      auto motherType = isK0S ? aod::femtodreamparticle::ParticleType::kV0K0Short
                              : aod::femtodreamparticle::ParticleType::kV0;
      outputPartsTime(timeStamp);
      outputPartsIndex(v0.globalIndex());
      outputParts(outputCollision.lastIndex(),
                  v0.pt(), v0.eta(), v0.phi(),
                  motherType,
                  cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kV0),
                  0,
                  v0.v0cosPA(),
                  indexChildID,
                  massV0,
                  antiMassV0);

      if (isDebug.value) {
        fillDebugParticle<true>(postrack);
        fillDebugParticle<true>(negtrack);
        fillDebugParticle<false>(v0, signV0);
      }

      if constexpr (IsMc) {
        fillMcParticle(col, v0, o2::aod::femtodreamparticle::ParticleType::kV0);
      }
      fIsV0Filled = true;
    }
    return fIsV0Filled;
  }

  template <DecayChannel Channel, bool IsMc, bool UseCharmMl, typename V0Type, typename CollisionType, typename TrackType, typename CandType>
  void fillCharmHadronTable(CollisionType const& col, TrackType const& tracks, V0Type const& fullV0s, CandType const& candidates)
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

    if (isNoSelectedV0s(col, tracks, fullV0s, k0sCuts) && sizeCand <= 0) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::RejNoV0sAndCharm);
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
    bool isV0Filled = false;
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
    isV0Filled = fillV0sForCharmHadron<IsMc>(col, fullV0s, tracks);
    aod::femtodreamcollision::BitMaskType bitV0 = 0;
    if (isV0Filled) {
      bitV0 |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::V0Selected);
    }

    aod::femtodreamcollision::BitMaskType bitCand = 0;
    if (sizeCand > 0) {
      bitCand |= 1 << 0;
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::CharmSelected);
    }

    if (isV0Filled && (sizeCand > 0)) {
      qaRegistry.fill(HIST("hEventQA"), 1 + Event::PairSelected);
    }

    rowMasks(bitV0,
             bitCand,
             0);
  }

  template <typename V0Type, typename TrackType>
  int determineV0Sign(V0Type const& v0, TrackType const& posTrack, TrackType const& negTrack)
  {
    const bool isK0S = (selectionFlagV0 == V0Channel::K0S);
    const bool isLambda = (selectionFlagV0 == V0Channel::Lambda);
    if (!isK0S && !isLambda) {
      LOG(fatal) << "Invalid V0 particle !! Please check the configuration";
    }
    // K0S self-conjugate: keep your convention
    if (isK0S) {
      return +1;
    }

    // Lambda / Anti-Lambda case

    const float mLamPDG = o2::constants::physics::MassLambda;
    const float mLamHyp = v0.mLambda();
    const float mAntiLamHyp = v0.mAntiLambda();

    const float diffLam = std::abs(mLamPDG - mLamHyp);
    const float diffAntiLam = std::abs(mLamPDG - mAntiLamHyp);

    const float offTPC = trkPIDnSigmaOffsetTPC.value;
    const float nSigmaPIDMax = trkPIDnSigmaForLambdaSign.value;

    // TPC n-sigma (apply offset by subtraction)
    const float prNeg = negTrack.tpcNSigmaPr() - offTPC;
    const float piPos = posTrack.tpcNSigmaPi() - offTPC;
    const float piNeg = negTrack.tpcNSigmaPi() - offTPC;
    const float prPos = posTrack.tpcNSigmaPr() - offTPC;

    const bool pidAntiLam = (std::abs(prNeg) < nSigmaPIDMax) && (std::abs(piPos) < nSigmaPIDMax);
    const bool pidLam = (std::abs(prPos) < nSigmaPIDMax) && (std::abs(piNeg) < nSigmaPIDMax);

    int sign = 0; // 0 = undecided

    // prefer: PID + mass preference (closer to PDG)
    if (pidAntiLam && (diffAntiLam < diffLam)) {
      sign = -1;
    } else if (pidLam && (diffLam <= diffAntiLam)) {
      sign = +1;
    } else {
      // fallback: PID only
      if (pidAntiLam) {
        sign = -1;
      } else {
        sign = +1;
      }
    }
    return sign;
  }

  template <typename TrackT>
  bool isV0DaughterPidSelected(const TrackT& posTrack,
                               const TrackT& negTrack,
                               int signV0,
                               bool& passPosPID,
                               bool& passNegPID)
  {
    passPosPID = false;
    passNegPID = false;

    const bool isK0S = (selectionFlagV0 == V0Channel::K0S);
    const bool isLambda = (selectionFlagV0 == V0Channel::Lambda);

    const float offTPC = trkPIDnSigmaOffsetTPC.value;

    const float prPos = posTrack.tpcNSigmaPr() - offTPC;
    const float piPos = posTrack.tpcNSigmaPi() - offTPC;
    const float prNeg = negTrack.tpcNSigmaPr() - offTPC;
    const float piNeg = negTrack.tpcNSigmaPi() - offTPC;

    if (isK0S) {
      const float nSigmaPiMaxK0S = V0Sel.confK0shortChildPIDnSigmaMax.value[0];
      passPosPID = (std::abs(piPos) < nSigmaPiMaxK0S);
      passNegPID = (std::abs(piNeg) < nSigmaPiMaxK0S);
      return passPosPID && passNegPID;
    }

    if (isLambda) {
      const float nSigmaPiMaxLam = V0Sel.confLambdaChildPIDnSigmaMax.value[0];
      const float nSigmaPrMaxLam = V0Sel.confLambdaChildPIDnSigmaMax.value[1];

      if (signV0 != +1 && signV0 != -1) {
        return false;
      }

      const bool posIsProton = (signV0 == +1);
      const bool negIsProton = (signV0 == -1);

      passPosPID = posIsProton ? (std::abs(prPos) < nSigmaPrMaxLam)
                               : (std::abs(piPos) < nSigmaPiMaxLam);

      passNegPID = negIsProton ? (std::abs(prNeg) < nSigmaPrMaxLam)
                               : (std::abs(piNeg) < nSigmaPiMaxLam);

      return passPosPID && passNegPID;
    }

    return false;
  }

  // check if there is no selected v0
  /// \param C type of the collision
  /// \param T type of the V0s
  /// \param TC type of the femto V0s cuts
  /// \return whether or not the V0s fulfills the all selections
  template <typename CollType, typename TrackType, typename TV0, typename CutType>
  bool isNoSelectedV0s(CollType const& col, TrackType const& /*col*/, TV0 const& V0s, CutType& v0Cuts)
  {
    for (auto const& v0 : V0s) {
      auto postrack = v0.template posTrack_as<TrackType>();
      auto negtrack = v0.template negTrack_as<TrackType>();
      if (v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
        return false;
      }
    }
    return true;
  }

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
                          o2::aod::V0Datas const& fullV0s,
                          soa::Filtered<CandidateD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, false, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataD0ToPiK, "Provide experimental data for D0ToPiK femto", false);

  void processDataD0ToPiKWithML(FemtoFullCollision const& col,
                                aod::BCsWithTimestamps const&,
                                FemtoHFTracks const& tracks,
                                o2::aod::V0Datas const& fullV0s,
                                soa::Filtered<soa::Join<CandidateD0,
                                                        aod::HfMlD0>> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, false, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataD0ToPiKWithML, "Provide experimental data for D0ToPiK with ml", false);

  void processMcD0ToPiK(FemtoFullCollisionMc const& col,
                        aod::BCsWithTimestamps const&,
                        FemtoHFMcTracks const& tracks,
                        soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                        aod::McParticles const&,
                        CandidateD0Mc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, true, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcD0ToPiK, "Provide Mc for D0ToPiK", false);

  void processMcD0ToPiKWithML(FemtoFullCollisionMc const& col,
                              aod::BCsWithTimestamps const&,
                              FemtoHFMcTracks const& tracks,
                              soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                              aod::McParticles const&,
                              soa::Join<CandidateD0Mc,
                                        aod::HfMlD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::D0ToPiK, true, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcD0ToPiKWithML, "Provide Mc for D0ToPiK with ml", false);

  void processMcD0ToPiKGen(Generated2ProngMc const& particles)
  {
    fillCharmHadMcGen<DecayChannel::D0ToPiK>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcD0ToPiKGen, "Provide Mc Generated D0ToPiK", false);

  /// DstarToD0Pi
  void processDataDstarToD0Pi(FemtoFullCollision const& col,
                              aod::BCsWithTimestamps const&,
                              FemtoHFTracks const& tracks,
                              o2::aod::V0Datas const& fullV0s,
                              soa::Filtered<CandidateDstar> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, false, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataDstarToD0Pi, "Provide experimental data for DstarToD0Pi femto", false);

  void processDataDstarToD0PiWithML(FemtoFullCollision const& col,
                                    aod::BCsWithTimestamps const&,
                                    FemtoHFTracks const& tracks,
                                    o2::aod::V0Datas const& fullV0s,
                                    soa::Filtered<soa::Join<CandidateDstar,
                                                            aod::HfMlDstarToD0Pi>> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, false, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataDstarToD0PiWithML, "Provide experimental data for DstarToD0Pi with ml", false);

  void processMcDstarToD0Pi(FemtoFullCollisionMc const& col,
                            aod::BCsWithTimestamps const&,
                            FemtoHFMcTracks const& tracks,
                            soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                            aod::McParticles const&,
                            CandidateDstarMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, true, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDstarToD0Pi, "Provide Mc for DstarToD0Pi", false);

  void processMcDstarToD0PiWithML(FemtoFullCollisionMc const& col,
                                  aod::BCsWithTimestamps const&,
                                  FemtoHFMcTracks const& tracks,
                                  soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                                  aod::McParticles const&,
                                  soa::Join<CandidateDstarMc,
                                            aod::HfMlDstarToD0Pi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DstarToD0Pi, true, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDstarToD0PiWithML, "Provide Mc for DstarToD0Pi with ml", false);

  void processMcDstarToD0PiGen(GeneratedDstarMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::DstarToD0Pi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDstarToD0PiGen, "Provide Mc Generated DstarToD0Pi", false);

  /// DplusToPiKPi
  void processDataDplusToPiKPi(FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               FemtoHFTracks const& tracks,
                               o2::aod::V0Datas const& fullV0s,
                               soa::Filtered<CandidateDplus> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, false, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataDplusToPiKPi, "Provide experimental data for DplusToPiKPi femto", false);

  void processDataDplusToPiKPiWithML(FemtoFullCollision const& col,
                                     aod::BCsWithTimestamps const&,
                                     FemtoHFTracks const& tracks,
                                     o2::aod::V0Datas const& fullV0s,
                                     soa::Filtered<soa::Join<CandidateDplus,
                                                             aod::HfMlDplusToPiKPi>> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, false, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataDplusToPiKPiWithML, "Provide experimental data for DplusToPiKPi with ml", false);

  void processMcDplusToPiKPi(FemtoFullCollisionMc const& col,
                             aod::BCsWithTimestamps const&,
                             FemtoHFMcTracks const& tracks,
                             soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                             aod::McParticles const&,
                             CandidateDplusMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, true, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDplusToPiKPi, "Provide Mc for DplusToPiKPi", false);

  void processMcDplusToPiKPiWithML(FemtoFullCollisionMc const& col,
                                   aod::BCsWithTimestamps const&,
                                   FemtoHFMcTracks const& tracks,
                                   soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                                   aod::McParticles const&,
                                   soa::Join<CandidateDplusMc,
                                             aod::HfMlDplusToPiKPi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::DplusToPiKPi, true, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDplusToPiKPiWithML, "Provide Mc for DplusToPiKPi with ml", false);

  void processMcDplusToPiKPiGen(Generated3ProngMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::DplusToPiKPi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcDplusToPiKPiGen, "Provide Mc Generated DplusToPiKPi", false);

  /// LcToPKPi
  void processDataLcToPKPi(FemtoFullCollision const& col,
                           aod::BCsWithTimestamps const&,
                           FemtoHFTracks const& tracks,
                           o2::aod::V0Datas const& fullV0s,
                           soa::Filtered<CandidateLc> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, false, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataLcToPKPi, "Provide experimental data for Lc(PKPi)-proton femto", false);

  void processDataLcToPKPiWithML(FemtoFullCollision const& col,
                                 aod::BCsWithTimestamps const&,
                                 FemtoHFTracks const& tracks,
                                 o2::aod::V0Datas const& fullV0s,
                                 soa::Filtered<soa::Join<CandidateLc,
                                                         aod::HfMlLcToPKPi>> const& candidates)
  {

    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, false, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processDataLcToPKPiWithML, "Provide experimental data for Lc(PKPi)-proton femto with ml", false);

  void processMcLcToPKPi(FemtoFullCollisionMc const& col,
                         aod::BCsWithTimestamps const&,
                         FemtoHFMcTracks const& tracks,
                         soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                         aod::McParticles const&,
                         CandidateLcMc const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, true, false>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcLcToPKPi, "Provide Mc for lctopkpi", false);

  void processMcLcToPKPiWithML(FemtoFullCollisionMc const& col,
                               aod::BCsWithTimestamps const&,
                               FemtoHFMcTracks const& tracks,
                               soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
                               aod::McParticles const&,
                               soa::Join<CandidateLcMc,
                                         aod::HfMlLcToPKPi> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    fillCharmHadronTable<DecayChannel::LcToPKPi, true, true>(col, tracks, fullV0s, candidates);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcLcToPKPiWithML, "Provide Mc for lctopkpi with ml", false);

  void processMcLcToPKPiGen(Generated3ProngMc const& particles)
  {

    fillCharmHadMcGen<DecayChannel::LcToPKPi>(particles);
  }
  PROCESS_SWITCH(HfProducerCharmHadronsV0FemtoDream, processMcLcToPKPiGen, "Provide Mc Generated lctopkpi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfProducerCharmHadronsV0FemtoDream>(cfgc)};
}
