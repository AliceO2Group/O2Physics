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

/// \file HFFilter.cxx
/// \brief task for selection of events with HF signals
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, INFN Bari

#include "EventFiltering/PWGHF/HFFilterHelpers.h"
#include "EventFiltering/filterTables.h"
//
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
//
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TString.h>

#include <Rtypes.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::analysis;
using namespace o2::aod::hffilters;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfFilter { // Main struct for HF triggers

  Produces<aod::HfFilters> tags;
  Produces<aod::HFOptimisationTreeBeauty> optimisationTreeBeauty;
  Produces<aod::HFOptimisationTreeCharm> optimisationTreeCharm;
  Produces<aod::HFOptimisationTreeFemto> optimisationTreeFemto;
  Produces<aod::HFOptimisationTreeCollisions> optimisationTreeCollisions;

  Configurable<int> activateQA{"activateQA", 0, "flag to enable QA histos (0 no QA, 1 basic QA, 2 extended QA, 3 very extended QA)"};
  Configurable<bool> activateSecVtxForB{"activateSecVtxForB", false, "flag to enable 2nd vertex fitting - only beauty hadrons"};

  // parameters for all triggers
  // event selection
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> applyTVX{"applyTVX", true, "flag to apply TVX selection"};
    Configurable<bool> applyTFBorderCut{"applyTFBorderCut", true, "apply TF border cut"};
    Configurable<bool> applyITSROFBorderCut{"applyITSROFBorderCut", false, "apply ITS ROF border cut"};
    Configurable<float> maxPvZ{"maxPvZ", 11.f, "maximum absolute value of PV position in the Z coordinate"};
  } evSel;

  // nsigma PID (except for V0 and cascades)
  Configurable<LabeledArray<float>> nSigmaPidCuts{"nSigmaPidCuts", {cutsNsigma[0], 4, 8, labelsRowsNsigma, labelsColumnsNsigma}, "Nsigma cuts for ITS/TPC/TOF PID (except for V0 and cascades)"};
  // min and max pts for tracks and bachelors (except for V0 and cascades)
  Configurable<LabeledArray<float>> ptCuts{"ptCuts", {cutsPt[0], 2, 10, labelsRowsCutsPt, labelsColumnsCutsPt}, "minimum and maximum pT for bachelor tracks (except for V0 and cascades)"};
  Configurable<LabeledArray<float>> trackQaulityCuts{"trackQaulityCuts", {cutsTrackQuality[0], 2, 7, labelsColumnsPtThresholdsForFemto, labelsColumnsTrackQuality}, "Track quality cuts for proton and deuteron)"};
  // parameters for high-pT triggers
  Configurable<LabeledArray<float>> ptThresholds{"ptThresholds", {cutsHighPtThresholds[0], 1, 2, labelsEmpty, labelsColumnsHighPtThresholds}, "pT treshold for high pT charm hadron candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCAXY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeautyToJPsi{"cutsTrackBeautyToJPsi", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for beauty->JPsi candidates (not muons)"};
  Configurable<std::string> paramCharmMassShape{"paramCharmMassShape", "2023_pass3", "Parametrisation of charm-hadron mass shape (options: 2023_pass3)"};
  Configurable<float> numSigmaDeltaMassCharmHad{"numSigmaDeltaMassCharmHad", 2.5, "Number of sigma for charm-hadron delta mass cut in B and D resonance triggers"};
  Configurable<std::vector<double>> pTBinsBHadron{"pTBinsBHadron", std::vector<double>{hf_trigger_cuts_presel_beauty::vecBinsPt}, "pT bin limits for beauty hadrons preselections"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<LabeledArray<double>> cutsBplus{"cutsBplus", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "B+ candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsBzeroToDstar{"cutsBzeroToDstar", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "B0 -> D*+ candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsBzero{"cutsBzero", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "B0 candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsBs{"cutsBs", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "Bs candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsBc{"cutsBc", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "Bc candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsLb{"cutsLb", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "Lb candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsXib{"cutsXib", {hf_trigger_cuts_presel_beauty::cuts[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVars, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsTopolBeauty}, "Xib candidate selection per pT bin"};
    Configurable<LabeledArray<double>> cutsBtoJPsiX{"cutsBtoJPsiX", {hf_trigger_cuts_presel_beauty::cutsBtoJPsi[0], hf_trigger_cuts_presel_beauty::nBinsPt, hf_trigger_cuts_presel_beauty::nCutVarsBtoJPsi, hf_trigger_cuts_presel_beauty::labelsPt, hf_trigger_cuts_presel_beauty::labelsColumnsCutsBeautyToJPsi}, "B->JPsiX candidate selection"};
  } cutsBtoHadrons;
  // parameters for femto triggers
  Configurable<float> femtoMaxRelativeMomentum{"femtoMaxRelativeMomentum", 2., "Maximal allowed value for relative momentum between charm-proton pairs in GeV/c"};
  Configurable<LabeledArray<int>> enableFemtoChannels{"enableFemtoChannels", {activeFemtoChannels[0], 2, 5, labelsRowsFemtoChannels, labelsColumnsFemtoChannels}, "Flags to enable/disable femto channels"};
  Configurable<LabeledArray<float>> ptThresholdsForFemto{"ptThresholdsForFemto", {cutsPtThresholdsForFemto[0], 1, 2, labelsEmpty, labelsColumnsPtThresholdsForFemto}, "pT treshold for proton or deuteron for kFemto triggers in GeV/c"};
  Configurable<bool> forceTofProtonForFemto{"forceTofProtonForFemto", true, "flag to force TOF PID for protons"};
  Configurable<bool> forceTofDeuteronForFemto{"forceTofDeuteronForFemto", false, "flag to force TOF PID for deuterons"};

  // double charm
  Configurable<LabeledArray<int>> enableDoubleCharmChannels{"enableDoubleCharmChannels", {activeDoubleCharmChannels[0], 2, 3, labelsRowsDoubleCharmChannels, labelsColumnsDoubleCharmChannels}, "Flags to enable/disable double charm channels"};
  Configurable<bool> keepOnlyDplusForDouble3Prongs{"keepOnlyDplusForDouble3Prongs", false, "Flag to enable/disable to keep only D+ in double charm 3-prongs trigger"};

  // parameters for resonance triggers
  Configurable<LabeledArray<float>> cutsGammaK0sLambda{"cutsGammaK0sLambda", {cutsV0s[0], 1, 6, labelsEmpty, labelsColumnsV0s}, "Selections for V0s (gamma, K0s, Lambda) for D+V0 triggers"};
  Configurable<LabeledArray<float>> cutsPtDeltaMassCharmReso{"cutsPtDeltaMassCharmReso", {cutsCharmReso[0], 4, 13, labelsRowsDeltaMassCharmReso, labelsColumnsDeltaMassCharmReso}, "pt (GeV/c) and invariant-mass delta (GeV/c2) for charm hadron resonances"};
  Configurable<bool> keepAlsoWrongDmesLambdaPairs{"keepAlsoWrongDmesLambdaPairs", true, "flat go keep also wrong sign D+Lambda pairs"};
  Configurable<bool> keepAlsoWrongDmesProtonPairs{"keepAlsoWrongDmesProtonPairs", true, "flat go keep also wrong sign D0p pairs"};
  Configurable<bool> keepAlsoWrongDstarMesProtonPairs{"keepAlsoWrongDstarMesProtonPairs", true, "flat go keep also wrong sign D*0p pairs"};

  // parameters for charm baryons to Xi bachelor
  Configurable<LabeledArray<float>> cutsXiCascades{"cutsXiCascades", {cutsCascades[0], 1, 8, labelsEmpty, labelsColumnsCascades}, "Selections for cascades (Xi) for Xi+bachelor triggers"};
  Configurable<LabeledArray<float>> cutsXiBachelor{"cutsXiBachelor", {cutsCharmBaryons[0], 1, 15, labelsEmpty, labelsColumnsCharmBarCuts}, "Selections for charm baryons (Xi+Pi, Xi+Ka, Xi+Pi+Pi)"};
  Configurable<LabeledArray<double>> cutsTrackCharmBaryonBachelor{"cutsTrackCharmBaryonBachelor", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for charm-baryon bachelor candidates"};
  Configurable<LabeledArray<int>> requireStrangenessTracking{"requireStrangenessTracking", {requireStrangenessTrackedXi[0], 1, 2, labelsEmpty, labelsColumnsCharmBaryons}, "Flags to require strangeness tracking for channels with Xi"};

  // parameters for ML application
  Configurable<std::vector<double>> pTBinsBDT{"pTBinsBDT", std::vector<double>{hf_cuts_bdt_multiclass::vecBinsPt}, "track pT bin limits for BDT cut"};

  Configurable<LabeledArray<double>> thresholdBDTScoreD0ToKPi{"thresholdBDTScoreD0ToKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D0 candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDPlusToPiKPi{"thresholdBDTScoreDPlusToPiKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of D+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreDSToPiKK{"thresholdBDTScoreDSToPiKK", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Ds+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreLcToPiKP{"thresholdBDTScoreLcToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreXicToPiKP{"thresholdBDTScoreXicToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Xic+ candidates"};

  Configurable<bool> acceptBdtBkgOnly{"acceptBdtBkgOnly", true, "Enable / disable selection based on BDT bkg score only"};

  // CCDB configuration
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // TPC PID calibrations
  Configurable<int> setTPCCalib{"setTPCCalib", 0, "0 is not use re-calibrations, 1 is compute TPC post-calibrated n-sigmas, 2 is using TPC Spline"};
  Configurable<std::string> ccdbBBProton{"ccdbBBProton", "Users/l/lserksny/PIDProton", "Path to the CCDB ocject for proton BB param"};
  Configurable<std::string> ccdbBBAntiProton{"ccdbBBAntiProton", "Users/l/lserksny/PIDAntiProton", "Path to the CCDB ocject for antiproton BB param"};
  Configurable<std::string> ccdbBBPion{"ccdbBBPion", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Pion BB param"};
  Configurable<std::string> ccdbBBAntiPion{"ccdbBBAntiPion", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiPion BB param"};
  Configurable<std::string> ccdbBBKaon{"ccdbBBKaon", "Users/l/lserksny/PIDPion", "Path to the CCDB ocject for Kaon BB param"};
  Configurable<std::string> ccdbBBAntiKaon{"ccdbBBAntiKaon", "Users/l/lserksny/PIDAntiPion", "Path to the CCDB ocject for antiKaon BB param"};
  Configurable<std::string> ccdbPathTPC{"ccdbPathTPC", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the CCDB object"};

  // parameter for Optimisation Tree
  Configurable<bool> applyOptimisation{"applyOptimisation", false, "Flag to enable or disable optimisation"};

  // manual downscale factors
  Configurable<bool> applyDownscale{"applyDownscale", false, "Flag to enable or disable the application of downscale factors"};
  Configurable<LabeledArray<double>> downscaleFactors{"downscaleFactors", {defDownscaleFactors[0], kNtriggersHF, 1, hfTriggerNames, labelsDownscaleFactor}, "Downscale factors for each trigger (from 0 to 1)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using BigTracksMCPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::McTrackLabels>;
  using BigTracksPID = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullDe, aod::pidTOFFullDe>;
  using TracksIUPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using CollsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;

  using Hf2ProngsWithMl = soa::Join<aod::Hf2Prongs, aod::Hf2ProngMlProbs>;
  using Hf3ProngsWithMl = soa::Join<aod::Hf3Prongs, aod::Hf3ProngMlProbs>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0s> v0sPerCollision = aod::v0::collisionId;
  Preslice<Hf2ProngsWithMl> hf2ProngPerCollision = aod::track_association::collisionId;
  Preslice<Hf3ProngsWithMl> hf3ProngPerCollision = aod::track_association::collisionId;
  Preslice<aod::Cascades> cascPerCollision = aod::cascade::collisionId;
  Preslice<aod::V0PhotonsKF> photonsPerCollision = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::AssignedTrackedCascades> trackedCascadesPerCollision = aod::track::collisionId;

  o2::ccdb::CcdbApi ccdbApi;
  int currentRun{0}; // needed to detect if the run changed and trigger update of calibrations etc.

  // array of BDT thresholds
  std::array<LabeledArray<double>, kNCharmParticles> thresholdBDTScores;

  o2::vertexing::DCAFitterN<2> df2;            // fitter for Charm Hadron vertex (2-prong vertex fitter)
  o2::vertexing::DCAFitterN<3> df3;            // fitter for Charm/Beauty Hadron vertex (3-prong vertex fitter)
  o2::vertexing::DCAFitterN<4> df4;            // fitter for Beauty Hadron vertex (4-prong vertex fitter)
  o2::vertexing::DCAFitterN<2> dfB;            // fitter for Beauty Hadron vertex (2-prong vertex fitter)
  o2::vertexing::DCAFitterN<3> dfBtoDstar;     // fitter for Beauty Hadron to D* vertex (3-prong vertex fitter)
  o2::vertexing::DCAFitterN<2> dfStrangeness;  // fitter for V0s and cascades (2-prong vertex fitter)
  o2::vertexing::DCAFitterN<3> dfStrangeness3; // fitter for Xic+ -> XiPiPi

  std::shared_ptr<TH1> hProcessedEvents;

  // QA histos
  std::shared_ptr<TH1> hN2ProngCharmCand, hN3ProngCharmCand;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmHighPt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmProtonKstarDistr{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmDeuteronKstarDistr{};
  std::array<std::shared_ptr<TH2>, nTotBeautyParts> hMassVsPtB{};
  std::array<std::shared_ptr<TH2>, kNCharmParticles + 23> hMassVsPtC{}; // +9 for resonances (D*+, D*0, Ds*+, Ds1+, Ds2*+, Xic+* right sign, Xic+* wrong sign, Xic0* right sign, Xic0* wrong sign) +2 for SigmaC (SigmaC++, SigmaC0) +2 for SigmaCK pairs (SigmaC++K-, SigmaC0K0s) +3 for charm baryons (Xi+Pi, Xi+Ka, Xi+Pi+Pi) + JPsi + 4 for charm baryons (D0+p, D0+pWrongSign, D*0p, D*0+pWrongSign)
  std::array<std::shared_ptr<TH2>, 4> hPrDePID;                         // proton TPC, proton TOF, deuteron TPC, deuteron TOF
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreBkg{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScorePrompt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreNonPrompt{};
  std::array<std::shared_ptr<TH2>, kNV0> hArmPod{};
  std::shared_ptr<TH2> hV0Selected;
  std::array<std::shared_ptr<TH2>, 2> hMassXi{}; // not tracked and tracked
  std::array<std::shared_ptr<TH2>, kNBeautyParticles> hCpaVsPtB{};
  std::array<std::shared_ptr<TH2>, kNBeautyParticles> hDecayLengthVsPtB{};
  std::array<std::shared_ptr<TH2>, kNBeautyParticles> hImpactParamProductVsPtB{};

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // helper object
  HfFilterHelper helper;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    helper.setHighPtTriggerThresholds(ptThresholds->get(0u, 0u), ptThresholds->get(0u, 1u));
    helper.setPtTriggerThresholdsForFemto(ptThresholdsForFemto->get(0u, 0u), ptThresholdsForFemto->get(0u, 1u));
    helper.setPtBinsSingleTracks(pTBinsTrack);
    helper.setPtBinsBeautyHadrons(pTBinsBHadron);
    helper.setPtLimitsBeautyBachelor(ptCuts->get(0u, 0u), ptCuts->get(1u, 0u), ptCuts->get(0u, 7u), ptCuts->get(1u, 7u));
    helper.setPtLimitsDstarSoftPion(ptCuts->get(0u, 1u), ptCuts->get(1u, 1u));
    helper.setPtLimitsProtonForFemto(ptCuts->get(0u, 2u), ptCuts->get(1u, 2u));
    helper.setPtLimitsDeuteronForFemto(ptCuts->get(0u, 6u), ptCuts->get(1u, 6u));
    helper.setPtLimitsCharmBaryonBachelor(ptCuts->get(0u, 3u), ptCuts->get(1u, 3u));
    helper.setPtLimitsLcResonanceBachelor(ptCuts->get(0u, 8u), ptCuts->get(1u, 8u));
    helper.setPtLimitsThetaCBachelor(ptCuts->get(0u, 9u), ptCuts->get(1u, 9u));
    helper.setCutsSingleTrackBeauty(cutsTrackBeauty3Prong, cutsTrackBeauty4Prong, cutsTrackBeauty4Prong);
    helper.setCutsSingleTrackCharmBaryonBachelor(cutsTrackCharmBaryonBachelor);
    helper.setCutsBhadrons(cutsBtoHadrons.cutsBplus, cutsBtoHadrons.cutsBzeroToDstar, cutsBtoHadrons.cutsBc, cutsBtoHadrons.cutsBzero, cutsBtoHadrons.cutsBs, cutsBtoHadrons.cutsLb, cutsBtoHadrons.cutsXib);
    helper.setCutsBtoJPsi(cutsBtoHadrons.cutsBtoJPsiX);
    helper.setNsigmaProtonCutsForFemto(std::array{nSigmaPidCuts->get(0u, 3u), nSigmaPidCuts->get(1u, 3u), nSigmaPidCuts->get(2u, 3u), nSigmaPidCuts->get(3u, 3u)});
    helper.setNsigmaDeuteronCutsForFemto(std::array{nSigmaPidCuts->get(0u, 6u), nSigmaPidCuts->get(1u, 6u), nSigmaPidCuts->get(2u, 6u), nSigmaPidCuts->get(3u, 6u)});
    helper.setDeuteronTrackSelectionForFemto(trackQaulityCuts->get(1u, 0u), trackQaulityCuts->get(1u, 1u), trackQaulityCuts->get(1u, 2u), trackQaulityCuts->get(1u, 3u), trackQaulityCuts->get(1u, 4u), trackQaulityCuts->get(1u, 5u), trackQaulityCuts->get(1u, 6u));
    helper.setNsigmaProtonCutsForCharmBaryons(nSigmaPidCuts->get(0u, 0u), nSigmaPidCuts->get(1u, 0u));
    helper.setNsigmaPionKaonCutsForDzero(nSigmaPidCuts->get(0u, 1u), nSigmaPidCuts->get(1u, 1u));
    helper.setNsigmaKaonCutsFor3Prongs(nSigmaPidCuts->get(0u, 2u), nSigmaPidCuts->get(1u, 2u));
    helper.setNsigmaKaonProtonCutsForBeautyToJPsi(nSigmaPidCuts->get(0u, 7u), nSigmaPidCuts->get(1u, 7u));
    helper.setForceTofForFemto(forceTofProtonForFemto, forceTofDeuteronForFemto);
    helper.setV0Selections(cutsGammaK0sLambda->get(0u, 0u), cutsGammaK0sLambda->get(0u, 1u), cutsGammaK0sLambda->get(0u, 2u), cutsGammaK0sLambda->get(0u, 3u), cutsGammaK0sLambda->get(0u, 4u), cutsGammaK0sLambda->get(0u, 5u));
    helper.setXiSelections(cutsXiCascades->get(0u, 0u), cutsXiCascades->get(0u, 1u), cutsXiCascades->get(0u, 2u), cutsXiCascades->get(0u, 3u), cutsXiCascades->get(0u, 4u), cutsXiCascades->get(0u, 5u), cutsXiCascades->get(0u, 6u), cutsXiCascades->get(0u, 7u));
    helper.setXiBachelorSelections(cutsXiBachelor->get(0u, 0u), cutsXiBachelor->get(0u, 1u), cutsXiBachelor->get(0u, 2u), cutsXiBachelor->get(0u, 3u), cutsXiBachelor->get(0u, 4u), cutsXiBachelor->get(0u, 5u), cutsXiBachelor->get(0u, 6u), cutsXiBachelor->get(0u, 7u), cutsXiBachelor->get(0u, 8u), cutsXiBachelor->get(0u, 9u), cutsXiBachelor->get(0u, 10u), cutsXiBachelor->get(0u, 11u), cutsXiBachelor->get(0u, 12u), cutsXiBachelor->get(0u, 13u), cutsXiBachelor->get(0u, 14u));
    helper.setNsigmaPiCutsForCharmBaryonBachelor(nSigmaPidCuts->get(0u, 4u), nSigmaPidCuts->get(1u, 4u));
    helper.setTpcPidCalibrationOption(setTPCCalib);
    helper.setMassResolParametrisation(paramCharmMassShape);
    helper.setNumSigmaForDeltaMassCharmHadCut(numSigmaDeltaMassCharmHad);
    helper.setPtRangeSoftPiSigmaC(ptCuts->get(0u, 4u), ptCuts->get(1u, 4u));
    helper.setPtDeltaMassRangeSigmaC(cutsPtDeltaMassCharmReso->get(0u, 6u), cutsPtDeltaMassCharmReso->get(1u, 6u), cutsPtDeltaMassCharmReso->get(0u, 7u), cutsPtDeltaMassCharmReso->get(1u, 7u), cutsPtDeltaMassCharmReso->get(0u, 8u), cutsPtDeltaMassCharmReso->get(1u, 8u), cutsPtDeltaMassCharmReso->get(0u, 9u), cutsPtDeltaMassCharmReso->get(1u, 9u), cutsPtDeltaMassCharmReso->get(2u, 6u), cutsPtDeltaMassCharmReso->get(2u, 7u), cutsPtDeltaMassCharmReso->get(2u, 8u), cutsPtDeltaMassCharmReso->get(2u, 9u));
    helper.setPtRangeSoftKaonXicResoToSigmaC(ptCuts->get(0u, 5u), ptCuts->get(1u, 5u));
    helper.setVtxConfiguration(dfStrangeness, true); // (DCAFitterN, useAbsDCA)
    helper.setVtxConfiguration(dfStrangeness3, true); // (DCAFitterN, useAbsDCA)
    dfStrangeness.setMatCorrType(matCorr);
    dfStrangeness3.setMatCorrType(matCorr);
    helper.setVtxConfiguration(df2, false); // (DCAFitterN, useAbsDCA)
    helper.setVtxConfiguration(df3, false);
    helper.setVtxConfiguration(df4, false);
    if (activateSecVtxForB) {
      helper.setVtxConfiguration(dfB, true);
      helper.setVtxConfiguration(dfBtoDstar, true);
    }

    // fetch config of track-index-skim-creator to apply the same cut on DeltaMassKK for Ds
    std::vector<double> ptBinsDsSkimCreator{};
    LabeledArray<double> cutsDsSkimCreator{};
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-track-index-skim-creator") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("binsPtDsToKKPi") == 0) {
            auto ptBins = option.defaultValue.get<double*>();
            double lastEl{-1.e6};
            int iPt{0};
            while (ptBins[iPt] > lastEl) {
              ptBinsDsSkimCreator.push_back(ptBins[iPt]);
              lastEl = ptBins[iPt];
              iPt++;
            }
          } else if (option.name.compare("cutsDsToKKPi") == 0) {
            cutsDsSkimCreator = option.defaultValue.get<LabeledArray<double>>();
          }
        }
      }
    }
    helper.setPreselDsToKKPi(ptBinsDsSkimCreator, cutsDsSkimCreator);

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;counts", HistType::kTH1D, {{kNtriggersHF + 2, -0.5, +kNtriggersHF + 1.5}});
    for (auto iBin = 0; iBin < kNtriggersHF + 2; ++iBin) {
      if (iBin < 2)
        hProcessedEvents->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        hProcessedEvents->GetXaxis()->SetBinLabel(iBin + 1, hfTriggerNames[iBin - 2].data());
    }

    if (activateQA) {
      hN2ProngCharmCand = registry.add<TH1>("fN2ProngCharmCand", "Number of 2-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1D, {{50, -0.5, 49.5}});
      hN3ProngCharmCand = registry.add<TH1>("fN3ProngCharmCand", "Number of 3-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1D, {{50, -0.5, 49.5}});
      for (int iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        hCharmHighPt[iCharmPart] = registry.add<TH1>(Form("f%sHighPt", charmParticleNames[iCharmPart].data()), Form("#it{p}_{T} distribution of triggered high-#it{p}_{T} %s candidates;#it{p}_{T} (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {ptAxis});
        hCharmProtonKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sProtonKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered p#minus%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {kstarAxis});
        hCharmDeuteronKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sDeuteronKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered de%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {kstarAxis});
        hMassVsPtC[iCharmPart] = registry.add<TH2>(Form("fMassVsPt%s", charmParticleNames[iCharmPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", charmParticleNames[iCharmPart].data()), HistType::kTH2D, {ptAxis, massAxisC[iCharmPart]});
        if (activateQA > 1) {
          hBDTScoreBkg[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreBkgDistr", charmParticleNames[iCharmPart].data()), Form("BDT background score distribution for %s;BDT background score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {bdtAxis});
          hBDTScorePrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScorePromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT prompt score distribution for %s;BDT prompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {bdtAxis});
          hBDTScoreNonPrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreNonPromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT nonprompt score distribution for %s;BDT nonprompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1D, {bdtAxis});
        }
      }
      // charm resonances
      hMassVsPtC[kNCharmParticles] = registry.add<TH2>("fMassVsPtDStarPlus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarPlus candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles]});
      hMassVsPtC[kNCharmParticles + 1] = registry.add<TH2>("fMassVsPtDStarZero", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarZero candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 1]});
      hMassVsPtC[kNCharmParticles + 2] = registry.add<TH2>("fMassVsPtDStarS", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered DStarS candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 2]});
      hMassVsPtC[kNCharmParticles + 3] = registry.add<TH2>("fMassVsPtDs1Plus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered Ds1Plus candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 3]});
      hMassVsPtC[kNCharmParticles + 4] = registry.add<TH2>("fMassVsPtDs2StarPlus", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered Ds2StarPlus candidates;#it{p}_{T} (GeV/#Delta#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 4]});
      hMassVsPtC[kNCharmParticles + 5] = registry.add<TH2>("fMassVsPtXicStarToDplusLambda", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered XicStar -> Dplus Lambda candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 5]});
      hMassVsPtC[kNCharmParticles + 6] = registry.add<TH2>("fMassVsPtXicStarToDplusLambdaWrongSign", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered opposite-sign XicStar -> Dplus Lambda candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 6]});
      hMassVsPtC[kNCharmParticles + 7] = registry.add<TH2>("fMassVsPtXicStarToD0Lambda", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered XicStar -> D0 Lambda candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 7]});
      hMassVsPtC[kNCharmParticles + 8] = registry.add<TH2>("fMassVsPtXicStarToD0LambdaWrongSign", "#Delta#it{M} vs. #it{p}_{T} distribution of triggered opposite-sign XicStar -> D0 Lambda candidates;#it{p}_{T} (GeV/#it{c});#Delta#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 8]});
      // SigmaC0,++
      hMassVsPtC[kNCharmParticles + 9] = registry.add<TH2>("fMassVsPtSigmaCPlusPlus", "#it{M}(pK#pi#pi)-M(pK#pi) vs. #it{p}_{T} distribution of #Sigma_{c}^{++} candidates for triggers;#it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});#it{M}(pK#pi#pi)-M(pK#pi);counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 9]});
      hMassVsPtC[kNCharmParticles + 10] = registry.add<TH2>("fMassVsPtSigmaC0", "#it{M}(pK#pi#pi)-M(pK#pi) vs. #it{p}_{T} distribution of #Sigma_{c}^{0} candidates for triggers;#it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});#it{M}(pK#pi#pi)-M(pK#pi);counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 10]});
      // SigmaCKaon pairs
      hMassVsPtC[kNCharmParticles + 11] = registry.add<TH2>("fMassVsPtSigmaC2455PlusPlusKaMinus", "#it{M}(#Sigma_{c}^{++}K^{-}(2455)) vs. #it{p}_{T} distribution of of triggered #Sigma_{c}^{++}K^{-} pairs;#it{p}_{T} (GeV/#it{c});#it{M}(#Sigma_{c}^{++}K^{-});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 11]});
      hMassVsPtC[kNCharmParticles + 12] = registry.add<TH2>("fMassVsPtSigmaC2520PlusPlusKaMinus", "#it{M}(#Sigma_{c}^{++}K^{-}(2520)) vs. #it{p}_{T} distribution of of triggered #Sigma_{c}^{++}K^{-} pairs;#it{p}_{T} (GeV/#it{c});#it{M}(#Sigma_{c}^{++}K^{-});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 12]});
      hMassVsPtC[kNCharmParticles + 13] = registry.add<TH2>("fMassVsPtSigmaC02455Ka0s", "#it{M}(#Sigma_{c}^{0}K^{0}_{s}(2455)) vs. #it{p}_{T} distribution of of triggered #Sigma_{c}^{0}K^{0}_{s} pairs;#it{p}_{T} (GeV/#it{c});#it{M}(#Sigma_{c}^{++}K^{-});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 13]});
      hMassVsPtC[kNCharmParticles + 14] = registry.add<TH2>("fMassVsPtSigmaC02520Ka0s", "#it{M}(#Sigma_{c}^{0}K^{0}_{s}(2520)) vs. #it{p}_{T} distribution of of triggered #Sigma_{c}^{0}K^{0}_{s} pairs;#it{p}_{T} (GeV/#it{c});#it{M}(#Sigma_{c}^{++}K^{-});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 14]});
      // charm baryons to LF cascades
      hMassVsPtC[kNCharmParticles + 15] = registry.add<TH2>("fMassVsPtCharmBaryonToXiPi", "#it{M} vs. #it{p}_{T} distribution of triggered #Xi+#pi candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 15]});
      hMassVsPtC[kNCharmParticles + 16] = registry.add<TH2>("fMassVsPtCharmBaryonToXiKa", "#it{M} vs. #it{p}_{T} distribution of triggered #Xi+K candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 16]});
      hMassVsPtC[kNCharmParticles + 17] = registry.add<TH2>("fMassVsPtCharmBaryonToXiPiPi", "#it{M} vs. #it{p}_{T} distribution of triggered #Xi+#pi+#pi candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 17]});
      // JPsi
      hMassVsPtC[kNCharmParticles + 18] = registry.add<TH2>("fMassVsPtJPsiToMuMu", "#it{M} vs. #it{p}_{T} distribution of triggered J/#psi to #mu#mu candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 18]});
      // Lc resonances
      hMassVsPtC[kNCharmParticles + 19] = registry.add<TH2>("fMassVsPtCharmBaryonToD0P", "#it{M} vs. #it{p}_{T} distribution of triggered D^{0}#p candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 19]});
      hMassVsPtC[kNCharmParticles + 20] = registry.add<TH2>("fMassVsPtCharmBaryonToD0PWrongSign", "#it{M} vs. #it{p}_{T} distribution of triggered D^{0}#p wrong sign candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 20]});
      // ThetaC
      hMassVsPtC[kNCharmParticles + 21] = registry.add<TH2>("fMassVsPtCharmBaryonToDstarP", "#it{M} vs. #it{p}_{T} distribution of triggered D^{*0}#p candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 21]});
      hMassVsPtC[kNCharmParticles + 22] = registry.add<TH2>("fMassVsPtCharmBaryonToDstarPWrongSign", "#it{M} vs. #it{p}_{T} distribution of triggered D^{*0}#p wrong sign candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {ptAxis, massAxisC[kNCharmParticles + 22]});

      for (int iBeautyPart{0}; iBeautyPart < kNBeautyParticles; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2D, {ptAxis, massAxisB[iBeautyPart]});
        hCpaVsPtB[iBeautyPart] = registry.add<TH2>(Form("fCpaVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("CPA vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2D, {ptAxis, {500, 0., 1}});
        hDecayLengthVsPtB[iBeautyPart] = registry.add<TH2>(Form("fDecayLengthVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("DecayLength vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2D, {ptAxis, {500, 0, 0.5}});
        if (iBeautyPart != kB0toDStar) {
          hImpactParamProductVsPtB[iBeautyPart] = registry.add<TH2>(Form("fImpactParamProductVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("ImpactParamProduct vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2D, {ptAxis, {500, -2.5e-3, +2.5e-3}});
        }
      }
      for (int iBeautyPart{kNBeautyParticles}; iBeautyPart < nTotBeautyParts; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2D, {ptAxis, massAxisB[iBeautyPart]});
      }
      constexpr int kNBinsHfVtxStages = kNHfVtxStage;
      std::string labels[kNBinsHfVtxStages];
      labels[HfVtxStage::Skimmed] = "Skimm CharmHad-Pi pairs";
      labels[HfVtxStage::BeautyVertex] = "vertex CharmHad-Pi pairs";
      labels[HfVtxStage::CharmHadPiSelected] = "selected CharmHad-Pi pairs";
      static const AxisSpec axisHfVtxStages = {kNBinsHfVtxStages, 0.5, kNBinsHfVtxStages + 0.5, ""};
      registry.add("fHfVtxStages", "HfVtxStages;;entries", HistType::kTH2D, {axisHfVtxStages, {kNBeautyParticles, -0.5, +kNBeautyParticles - 0.5}});
      for (int iBin = 0; iBin < kNBinsHfVtxStages; iBin++) {
        registry.get<TH2>(HIST("fHfVtxStages"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
      for (int iBin = 0; iBin < kNBeautyParticles; iBin++) {
        registry.get<TH2>(HIST("fHfVtxStages"))->GetYaxis()->SetBinLabel(iBin + 1, beautyParticleNames[iBin].data());
      }

      for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
        hArmPod[iV0] = registry.add<TH2>(Form("fArmPod%s", v0Names[iV0].data()), Form("Armenteros Podolanski plot for selected %s;#it{#alpha};#it{q}_{T} (GeV/#it{c})", v0Labels[iV0].data()), HistType::kTH2D, {alphaAxis, qtAxis});
      }
      hMassXi[0] = registry.add<TH2>("fMassXi", "#it{M} distribution of #Xi candidates;sign;#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {{3, -1.5f, 1.5f}, {100, 1.28f, 1.36f}});
      hMassXi[1] = registry.add<TH2>("fMassTrackedXi", "#it{M} distribution of #Xi candidates;sign;#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2D, {{3, -1.5f, 1.5f}, {100, 1.28f, 1.36f}});

      if (activateQA > 1) {
        hPrDePID[0] = registry.add<TH2>("fProtonTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2D, {pAxis, nSigmaAxis});
        hPrDePID[1] = registry.add<TH2>("fProtonTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2D, {pAxis, nSigmaAxis});
        hPrDePID[2] = registry.add<TH2>("fDeuteronTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected deuterons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2D, {pAxis, nSigmaAxis});
        hPrDePID[3] = registry.add<TH2>("fDeuteronTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected deuterons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2D, {pAxis, nSigmaAxis});

        hV0Selected = registry.add<TH2>("fV0Selected", "Selections for V0s;;counts", HistType::kTH2D, {{9, -0.5, 8.5}, {kNV0, -0.5, +kNV0 - 0.5}});

        for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
          hV0Selected->GetYaxis()->SetBinLabel(iV0 + 1, v0Labels[iV0].data());
        }
        hV0Selected->GetXaxis()->SetBinLabel(1, "analysed");
        hV0Selected->GetXaxis()->SetBinLabel(2, "rej. |#eta|");
        hV0Selected->GetXaxis()->SetBinLabel(3, "rej. radius");
        hV0Selected->GetXaxis()->SetBinLabel(4, "rej. cos(#theta_{P})");
        hV0Selected->GetXaxis()->SetBinLabel(5, "rej. Mass");
        hV0Selected->GetXaxis()->SetBinLabel(6, "rej. DCA V0");
        hV0Selected->GetXaxis()->SetBinLabel(7, "rej. DCA V0 daughters");
        hV0Selected->GetXaxis()->SetBinLabel(8, "rej. PID");
        hV0Selected->GetXaxis()->SetBinLabel(9, "selected");
      }
    }

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    thresholdBDTScores = {thresholdBDTScoreD0ToKPi, thresholdBDTScoreDPlusToPiKPi, thresholdBDTScoreDSToPiKK, thresholdBDTScoreLcToPiKP, thresholdBDTScoreXicToPiKP};
  }

  void process(CollsWithEvSel const& collisions,
               aod::BCsWithTimestamps const&,
               aod::V0s const& v0s,
               aod::Cascades const& cascades,
               aod::AssignedTrackedCascades const& trackedCasc,
               Hf2ProngsWithMl const& cand2Prongs,
               Hf3ProngsWithMl const& cand3Prongs,
               aod::TrackAssoc const& trackIndices,
               BigTracksPID const& tracks,
               TracksIUPID const& tracksIU,
               aod::V0PhotonsKF const& photons,
               aod::V0Legs const&)
  {
    for (const auto& collision : collisions) {

      bool keepEvent[kNtriggersHF]{false};
      bool isSelectedTVX = evSel.applyTVX ? collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) : true;
      bool isSelectedTFBorder = evSel.applyTFBorderCut ? collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) : true;
      bool isSelectedITSROFBorder = evSel.applyITSROFBorderCut ? collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) : true;
      bool isSelectedPvZ = (std::fabs(collision.posZ()) < evSel.maxPvZ);
      if (!isSelectedTVX || !isSelectedTFBorder || !isSelectedITSROFBorder || !isSelectedPvZ) {
        tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix], keepEvent[kV0Charm2P], keepEvent[kV0Charm3P], keepEvent[kCharmBarToXiBach], keepEvent[kSigmaCPPK], keepEvent[kSigmaC0K0], keepEvent[kPhotonCharm2P], keepEvent[kPhotonCharm3P], keepEvent[kSingleCharm2P], keepEvent[kSingleCharm3P], keepEvent[kSingleNonPromptCharm2P], keepEvent[kSingleNonPromptCharm3P], keepEvent[kCharmBarToXi2Bach], keepEvent[kPrCharm2P], keepEvent[kBtoJPsiKa], keepEvent[kBtoJPsiKstar], keepEvent[kBtoJPsiPhi], keepEvent[kBtoJPsiPrKa], keepEvent[kBtoJPsiPi]);
        continue;
      }

      auto thisCollId = collision.globalIndex();

      if (applyOptimisation) {
        optimisationTreeCollisions(thisCollId);
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      // needed for track propagation
      if (currentRun != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        // setMatLUT only after magfield has been initalized
        // (setMatLUT has implicit and problematic init field call if not)
        o2::base::Propagator::Instance()->setMatLUT(lut);

        // needed for TPC PID postcalibrations
        if (setTPCCalib == 1) {
          helper.setTpcRecalibMaps(ccdb, bc, ccdbPathTPC);
        } else if (setTPCCalib > 1) {
          helper.setValuesBB(ccdbApi, bc, std::array{ccdbBBPion.value, ccdbBBAntiPion.value, ccdbBBKaon.value, ccdbBBAntiKaon.value, ccdbBBProton.value, ccdbBBAntiProton.value, ccdbBBProton.value, ccdbBBAntiProton.value}); // dummy for deuteron
        }

        auto bz = o2::base::Propagator::Instance()->getNominalBz();
        dfStrangeness.setBz(bz);
        dfStrangeness3.setBz(bz);
        df2.setBz(bz);
        df3.setBz(bz);
        df4.setBz(bz);
        if (activateSecVtxForB) {
          dfB.setBz(bz);
          dfBtoDstar.setBz(bz);
        }

        currentRun = bc.runNumber();
      }

      hProcessedEvents->Fill(0);

      std::vector<std::vector<int64_t>> indicesDau2Prong{}, indicesDau2ProngPrompt{};

      auto cand2ProngsThisColl = cand2Prongs.sliceBy(hf2ProngPerCollision, thisCollId);
      for (const auto& cand2Prong : cand2ProngsThisColl) { // start loop over 2 prongs

        int8_t preselD0 = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK);            // check if it's a D0
        int8_t preselJPsiToMuMu = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::JpsiToMuMu); // check if it's a JPsi
        if (preselD0 == 0 && preselJPsiToMuMu == 0) {
          continue;
        }

        auto trackPos = tracks.rawIteratorAt(cand2Prong.prong0Id()); // positive daughter
        auto trackNeg = tracks.rawIteratorAt(cand2Prong.prong1Id()); // negative daughter

        if (preselD0) {
          preselD0 = helper.isDzeroPreselected(trackPos, trackNeg);
        }

        auto trackParPos = getTrackParCov(trackPos);
        auto trackParNeg = getTrackParCov(trackNeg);
        std::array<float, 2> dcaPos{trackPos.dcaXY(), trackPos.dcaZ()};
        std::array<float, 2> dcaNeg{trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::array<float, 3> pVecPos{trackPos.pVector()};
        std::array<float, 3> pVecNeg{trackNeg.pVector()};
        if (trackPos.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPos, 2.f, noMatCorr, &dcaPos);
          getPxPyPz(trackParPos, pVecPos);
        }
        if (trackNeg.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParNeg, 2.f, noMatCorr, &dcaNeg);
          getPxPyPz(trackParNeg, pVecNeg);
        }

        // apply ML models for D0
        bool isD0CharmTagged{false}, isD0BeautyTagged{false}, isD0SignalTagged{false};
        std::vector<float> scores{};
        if (preselD0) {
          scores.insert(scores.end(), cand2Prong.mlProbSkimD0ToKPi().begin(), cand2Prong.mlProbSkimD0ToKPi().end());
          if (scores.size() != 3) {
            scores.resize(3);
            scores[0] = 2.;
            scores[1] = -1.;
            scores[2] = -1.;
          }
          auto tagBDT = helper.isBDTSelected(scores, thresholdBDTScores[kD0]);
          isD0CharmTagged = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
          isD0BeautyTagged = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);
          isD0SignalTagged = acceptBdtBkgOnly ? TESTBIT(tagBDT, RecoDecay::OriginType::None) : (isD0CharmTagged || isD0BeautyTagged);

          if (activateQA > 1) {
            hBDTScoreBkg[kD0]->Fill(scores[0]);
            hBDTScorePrompt[kD0]->Fill(scores[1]);
            hBDTScoreNonPrompt[kD0]->Fill(scores[2]);
          }
        }

        auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
        auto pt2Prong = RecoDecay::pt(pVec2Prong);

        if (preselJPsiToMuMu) {
          float ptMuonMin = cutsBtoHadrons.cutsBtoJPsiX->get(0u, 0u); // assuming that the cut is looser in the first pT bin
          auto ptPos = RecoDecay::pt(pVecPos);
          auto ptNeg = RecoDecay::pt(pVecNeg);
          if (ptPos < ptMuonMin || ptNeg < ptMuonMin) {
            preselJPsiToMuMu = 0u;
          } else {
            auto massJPsiCand = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massMu, massMu});
            hMassVsPtC[kNCharmParticles + 18]->Fill(pt2Prong, massJPsiCand);
          }
        }

        if (!isD0SignalTagged && !preselJPsiToMuMu) {
          continue;
        }

        int8_t selD0InMass{0};
        double massD0Cand{-1.}, massD0BarCand{-1.};
        if (isD0SignalTagged) {
          // single D0
          keepEvent[kSingleCharm2P] = true;
          if (isD0BeautyTagged) {
            keepEvent[kSingleNonPromptCharm2P] = true;
          }
          // single D0 at high pT
          if (helper.isSelectedHighPt2Prong(pt2Prong)) {
            keepEvent[kHighPt2P] = true;
            if (activateQA) {
              hCharmHighPt[kD0]->Fill(pt2Prong);
            }
          }
          // multi-charm selection
          indicesDau2Prong.push_back(std::vector<int64_t>{trackPos.globalIndex(), trackNeg.globalIndex()});
          if (isD0CharmTagged) {
            indicesDau2ProngPrompt.push_back(std::vector<int64_t>{trackPos.globalIndex(), trackNeg.globalIndex()});
          }

          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, o2::constants::physics::Pdg::kD0, pt2Prong, scores[0], scores[1], scores[2]);
          }
          selD0InMass = helper.isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong, preselD0, activateQA, hMassVsPtC[kD0]);
          // compute masses already here, needed both for B0 --> D* (--> D0 Pi) Pi and Ds1 --> D* (--> D0 Pi) K0S
          massD0Cand = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massKa});
          massD0BarCand = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massKa, massPi});
        }

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        auto tracksWithItsPid = soa::Attach<BigTracksPID, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
        for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks
          auto track = tracksWithItsPid.rawIteratorAt(trackId.trackId());

          if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
            continue;
          }

          auto trackParThird = getTrackParCov(track);
          std::array<float, 2> dcaThird{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecThird = track.pVector();
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
            getPxPyPz(trackParThird, pVecThird);
          }

          // Beauty with D0
          if (!keepEvent[kBeauty3P] && isD0BeautyTagged) {
            int16_t isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty<kBeauty3P>(track, trackParThird, dcaThird);
            if (TESTBIT(isTrackSelected, kForBeauty) && ((TESTBIT(selD0InMass, 0) && track.sign() < 0) || (TESTBIT(selD0InMass, 1) && track.sign() > 0))) { // D0 pi-/K- and D0bar pi+/K+
              auto massCandD0Pi = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
              auto massCandD0K = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massKa});
              auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
              auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
              bool isBplusInMass = helper.isSelectedBhadronInMassRange(ptCand, massCandD0Pi, kBplus);
              bool isBcInMass = helper.isSelectedBhadronInMassRange(ptCand, massCandD0K, kBc);

              if (TESTBIT(isTrackSelected, kForBeauty) && (isBplusInMass || isBcInMass)) {
                if (activateQA) {
                  if (isBplusInMass)
                    registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::Skimmed, kBplus);
                  if (isBcInMass)
                    registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::Skimmed, kBc);
                }
                if (!activateSecVtxForB) {
                  keepEvent[kBeauty3P] = true;
                  // fill optimisation tree for D0
                  if (applyOptimisation) {
                    optimisationTreeBeauty(thisCollId, o2::constants::physics::Pdg::kD0, pt2Prong, scores[0], scores[1], scores[2], dcaThird[0]);
                  }
                  if (activateQA) {
                    if (isBplusInMass)
                      hMassVsPtB[kBplus]->Fill(ptCand, massCandD0Pi);
                    if (isBcInMass)
                      hMassVsPtB[kBc]->Fill(ptCand, massCandD0K);
                  }
                } else {
                  int nVtxD{0};
                  try {
                    nVtxD = df2.process(trackParPos, trackParNeg);
                  } catch (...) {
                    LOG(error) << "Exception caught in DCA fitter process call for charm 2-prong!";
                    nVtxD = 0;
                  }
                  if (nVtxD != 0) {
                    std::array<float, 3> pVecPosVtx{}, pVecNegVtx{};
                    df2.getTrack(0).getPxPyPzGlo(pVecPosVtx);
                    df2.getTrack(1).getPxPyPzGlo(pVecNegVtx);
                    auto trackParD = df2.createParentTrackParCov();
                    trackParD.setAbsCharge(0); // to be sure
                    auto pVec2ProngVtx = RecoDecay::pVec(pVecPosVtx, pVecNegVtx);
                    int nVtxB{0};
                    try {
                      nVtxB = dfB.process(trackParD, trackParThird);
                    } catch (...) {
                      LOG(error) << "Exception caught in DCA fitter process call for beauty 3-prong!";
                      nVtxB = 0;
                    }
                    if (nVtxB != 0) {
                      if (activateQA) {
                        registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::BeautyVertex, kBplus);
                      }
                      const auto& secondaryVertexBtoD0h = dfB.getPCACandidate();
                      std::array<float, 3> pVecThirdVtx{};
                      dfB.getTrack(0).getPxPyPzGlo(pVec2ProngVtx);
                      dfB.getTrack(1).getPxPyPzGlo(pVecThirdVtx);
                      std::array<float, 2> dca2Prong; //{trackParD.dcaXY(), trackParD.dcaZ()};
                      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParD, 2.f, noMatCorr, &dca2Prong);
                      bool isBplus = helper.isSelectedBhadron(pVec2ProngVtx, pVecThirdVtx, dca2Prong, dcaThird, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}, kBplus);
                      bool isBc = helper.isSelectedBhadron(pVec2ProngVtx, pVecThirdVtx, dca2Prong, dcaThird, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}, kBc);

                      if (isBplus || isBc) {
                        keepEvent[kBeauty3P] = true;
                        // fill optimisation tree for D0
                        if (applyOptimisation) {
                          optimisationTreeBeauty(thisCollId, o2::constants::physics::Pdg::kD0, pt2Prong, scores[0], scores[1], scores[2], dcaThird[0]);
                        }
                        if (activateQA) {
                          if (isBplus) {
                            registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::CharmHadPiSelected, kBplus);
                            hCpaVsPtB[kBplus]->Fill(ptCand, RecoDecay::cpa(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}, RecoDecay::pVec(pVec2ProngVtx, pVecThirdVtx)));
                            hDecayLengthVsPtB[kBplus]->Fill(ptCand, RecoDecay::distance(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}));
                            hImpactParamProductVsPtB[kBplus]->Fill(ptCand, dca2Prong[0] * dcaThird[0]);
                            hMassVsPtB[kBplus]->Fill(ptCand, massCandD0Pi);
                          }
                          if (isBc) {
                            registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::CharmHadPiSelected, kBc);
                            hCpaVsPtB[kBc]->Fill(ptCand, RecoDecay::cpa(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}, RecoDecay::pVec(pVec2ProngVtx, pVecThirdVtx)));
                            hDecayLengthVsPtB[kBc]->Fill(ptCand, RecoDecay::distance(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBtoD0h[0], secondaryVertexBtoD0h[1], secondaryVertexBtoD0h[2]}));
                            hImpactParamProductVsPtB[kBc]->Fill(ptCand, dca2Prong[0] * dcaThird[0]);
                            hMassVsPtB[kBc]->Fill(ptCand, massCandD0K);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            if (!keepEvent[kBeauty3P] && TESTBIT(isTrackSelected, kSoftPionForBeauty) && ((TESTBIT(selD0InMass, 0) && track.sign() > 0) || (TESTBIT(selD0InMass, 1) && track.sign() < 0))) { // D0 pi+ and D0bar pi-
              auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
              auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
              std::array<float, 2> massDausD0{massPi, massKa};
              auto massD0dau = massD0Cand;
              if (track.sign() < 0) {
                massDausD0[0] = massKa;
                massDausD0[1] = massPi;
                massD0dau = massD0BarCand;
              }
              auto massDstarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecThird}, std::array{massDausD0[0], massDausD0[1], massPi});
              auto massDiffDstar = massDstarCand - massD0dau;
              if (cutsPtDeltaMassCharmReso->get(0u, 0u) <= massDiffDstar && massDiffDstar <= cutsPtDeltaMassCharmReso->get(1u, 0u) && ptCand > cutsPtDeltaMassCharmReso->get(2u, 0u)) { // additional check for B0->D*pi polarization studies
                if (activateQA) {
                  hMassVsPtC[kNCharmParticles]->Fill(ptCand, massDiffDstar);
                }
                for (const auto& trackIdB : trackIdsThisCollision) { // start loop over tracks
                  auto trackB = tracks.rawIteratorAt(trackIdB.trackId());
                  if (track.globalIndex() == trackB.globalIndex()) {
                    continue;
                  }
                  auto trackParFourth = getTrackParCov(trackB);
                  std::array<float, 2> dcaFourth{trackB.dcaXY(), trackB.dcaZ()};
                  std::array<float, 3> pVecFourth = trackB.pVector();
                  if (trackB.collisionId() != thisCollId) {
                    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFourth, 2.f, noMatCorr, &dcaFourth);
                    getPxPyPz(trackParFourth, pVecFourth);
                  }

                  auto isTrackFourthSelected = helper.isSelectedTrackForSoftPionOrBeauty<kBeauty3P>(trackB, trackParFourth, dcaFourth);
                  if (track.sign() * trackB.sign() < 0 && TESTBIT(isTrackFourthSelected, kForBeauty)) {
                    auto massCandB0 = RecoDecay::m(std::array{pVecBeauty3Prong, pVecFourth}, std::array{massDStar, massPi});
                    auto pVecBeauty4Prong = RecoDecay::pVec(pVec2Prong, pVecThird, pVecFourth);
                    auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                    if (helper.isSelectedBhadronInMassRange(ptCandBeauty4Prong, massCandB0, kB0toDStar)) {
                      if (activateQA) {
                        registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::Skimmed, kB0toDStar);
                      }
                      if (!activateSecVtxForB) {
                        keepEvent[kBeauty3P] = true;
                        // fill optimisation tree for D*
                        if (applyOptimisation) {
                          optimisationTreeBeauty(thisCollId, 413, pt2Prong, scores[0], scores[1], scores[2], dcaFourth[0]); // pdgCode of D*(2010)+: 413
                        }
                        if (activateQA) {
                          hMassVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, massCandB0);
                        }
                      } else {
                        int nVtxD{0};
                        try {
                          nVtxD = df2.process(trackParPos, trackParNeg);
                        } catch (...) {
                          LOG(error) << "Exception caught in DCA fitter process call for charm 2-prong!";
                          nVtxD = 0;
                        }
                        if (nVtxD > 0) {
                          std::array<float, 3> pVecPosVtx{}, pVecNegVtx{};
                          df2.getTrack(0).getPxPyPzGlo(pVecPosVtx);
                          df2.getTrack(1).getPxPyPzGlo(pVecNegVtx);
                          auto trackParD = df2.createParentTrackParCov();
                          trackParD.setAbsCharge(0); // to be sure
                          auto pVec2ProngVtx = RecoDecay::pVec(pVecPosVtx, pVecNegVtx);
                          int nVtxB{0};
                          try {
                            nVtxB = dfBtoDstar.process(trackParD, trackParThird, trackParFourth);
                          } catch (...) {
                            LOG(error) << "Exception caught in DCA fitter process call for beauty to D*!";
                            nVtxB = 0;
                          }
                          if (nVtxB != 0) {
                            if (activateQA) {
                              registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::BeautyVertex, kB0toDStar);
                            }
                            const auto& secondaryVertexBzero = dfBtoDstar.getPCACandidate();
                            std::array<float, 3> pVecThirdVtx{}, pVecFourthVtx{};
                            dfBtoDstar.getTrack(0).getPxPyPzGlo(pVec2ProngVtx);
                            dfBtoDstar.getTrack(1).getPxPyPzGlo(pVecThirdVtx);
                            dfBtoDstar.getTrack(2).getPxPyPzGlo(pVecFourthVtx);
                            bool isBzero = helper.isSelectedBzeroToDstar(pVec2ProngVtx, pVecThirdVtx, pVecFourthVtx, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBzero[0], secondaryVertexBzero[1], secondaryVertexBzero[2]});
                            if (isBzero) {
                              keepEvent[kBeauty3P] = true;
                              // fill optimisation tree for D0
                              if (applyOptimisation) {
                                optimisationTreeBeauty(thisCollId, 413, pt2Prong, scores[0], scores[1], scores[2], dcaFourth[0]); // pdgCode of D*(2010)+: 413
                              }
                              if (activateQA) {
                                registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::CharmHadPiSelected, kB0toDStar);
                                hCpaVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, RecoDecay::cpa(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBzero[0], secondaryVertexBzero[1], secondaryVertexBzero[2]}, RecoDecay::pVec(pVec2ProngVtx, pVecThirdVtx, pVecFourthVtx)));
                                hDecayLengthVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, RecoDecay::distance(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBzero[0], secondaryVertexBzero[1], secondaryVertexBzero[2]}));
                                hMassVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, massCandB0);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } // end beauty selection

            // 2-prong femto
            bool isProtonForCharm2Prong = helper.isSelectedTrack4Femto(track, trackParThird, activateQA, hPrDePID[0], hPrDePID[1], kProtonForFemto);
            bool isDeuteronForCharm2Prong = helper.isSelectedTrack4Femto(track, trackParThird, activateQA, hPrDePID[2], hPrDePID[3], kDeuteronForFemto);

            if (track.collisionId() == thisCollId) {
              if (isProtonForCharm2Prong && !keepEvent[kFemto2P] && enableFemtoChannels->get(0u, 0u) && isD0CharmTagged) {
                float relativeMomentum = helper.computeRelativeMomentum(pVecThird, pVec2Prong, massD0);
                if (applyOptimisation) {
                  optimisationTreeFemto(thisCollId, o2::constants::physics::Pdg::kD0, pt2Prong, scores[0], scores[1], scores[2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr(), track.tpcNSigmaDe(), track.tofNSigmaDe());
                }
                if (relativeMomentum < femtoMaxRelativeMomentum) {
                  keepEvent[kFemto2P] = true;
                  if (activateQA) {
                    hCharmProtonKstarDistr[kD0]->Fill(relativeMomentum);
                  }
                }
              }

              if (isDeuteronForCharm2Prong && !keepEvent[kFemto2P] && enableFemtoChannels->get(1u, 0u) && isD0CharmTagged) {
                float relativeMomentum = helper.computeRelativeMomentum(pVecThird, pVec2Prong, massD0);
                if (applyOptimisation) {
                  optimisationTreeFemto(thisCollId, o2::constants::physics::Pdg::kD0, pt2Prong, scores[0], scores[1], scores[2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr(), track.tpcNSigmaDe(), track.tofNSigmaDe());
                }
                if (relativeMomentum < femtoMaxRelativeMomentum) {
                  keepEvent[kFemto2P] = true;
                  if (activateQA) {
                    hCharmDeuteronKstarDistr[kD0]->Fill(relativeMomentum);
                  }
                }
              }
            } // end femto charm 2prong selection

          // Beauty with JPsi
          if (preselJPsiToMuMu) {
            if (!TESTBIT(helper.isSelectedTrackForSoftPionOrBeauty<kBtoJPsiKa>(track, trackParThird, dcaThird), kForBeauty)) { // same for all channels
              continue;
            }
            std::array<float, 3> pVecPosVtx{}, pVecNegVtx{}, pVecThirdVtx{}, pVecFourthVtx{};
            // 3-prong vertices
            if (!keepEvent[kBtoJPsiKa] || !keepEvent[kBtoJPsiPi]) {
              int nVtxB{0};
              try {
                nVtxB = df3.process(trackParPos, trackParNeg, trackParThird);
              } catch (...) {
                LOG(error) << "Exception caught in DCA fitter process call for beauty to JPsi 3-prong!";
                nVtxB = 0;
              }
              if (nVtxB != 0) {
                const auto& secondaryVertexBto3tracks = df3.getPCACandidate();
                df3.getTrack(0).getPxPyPzGlo(pVecPosVtx);
                df3.getTrack(1).getPxPyPzGlo(pVecNegVtx);
                df3.getTrack(2).getPxPyPzGlo(pVecThirdVtx);
                auto isBhadSel = helper.isSelectedBhadronToJPsi<3>(std::array{pVecPosVtx, pVecNegVtx, pVecThirdVtx}, std::array{track}, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBto3tracks[0], secondaryVertexBto3tracks[1], secondaryVertexBto3tracks[2]}, activateQA, hMassVsPtB);
                if (TESTBIT(isBhadSel, kBplusToJPsi)) {
                  keepEvent[kBtoJPsiKa] = true;
                }
                if (TESTBIT(isBhadSel, kBcToJPsi)) {
                  keepEvent[kBtoJPsiPi] = true;
                }
              }
            }
            // 4-prong vertices
            if (!keepEvent[kBtoJPsiKstar] || !keepEvent[kBtoJPsiPhi] || !keepEvent[kBtoJPsiPrKa]) {
              for (const auto& trackIdB : trackIdsThisCollision) { // start loop over tracks
                if (keepEvent[kBtoJPsiKstar] && keepEvent[kBtoJPsiPhi] && keepEvent[kBtoJPsiPrKa]) {
                  break;
                }
                auto trackFourth = tracksWithItsPid.rawIteratorAt(trackIdB.trackId());
                if (trackFourth.globalIndex() == track.globalIndex() || trackFourth.globalIndex() == trackPos.globalIndex() || trackFourth.globalIndex() == trackNeg.globalIndex() || trackFourth.sign() * track.sign() > 0) {
                  continue;
                }
                auto trackParFourth = getTrackParCov(trackFourth);
                std::array<float, 2> dcaFourth{trackFourth.dcaXY(), trackFourth.dcaZ()};
                std::array<float, 3> pVecFourth = trackFourth.pVector();
                if (trackFourth.collisionId() != thisCollId) {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFourth, 2.f, noMatCorr, &dcaFourth);
                  getPxPyPz(trackParFourth, pVecFourth);
                }
                if (!TESTBIT(helper.isSelectedTrackForSoftPionOrBeauty<kBtoJPsiKa>(trackFourth, trackParFourth, dcaFourth), kForBeauty)) { // same for all channels
                  continue;
                }
                int nVtxB{0};
                try {
                  nVtxB = df4.process(trackParPos, trackParNeg, trackParThird, trackParFourth);
                } catch (...) {
                  LOG(error) << "Exception caught in DCA fitter process call for beauty to JPsi 4-prong!";
                  nVtxB = 0;
                }
                if (nVtxB != 0) {
                  const auto& secondaryVertexBto4tracks = df4.getPCACandidate();
                  df4.getTrack(0).getPxPyPzGlo(pVecPosVtx);
                  df4.getTrack(1).getPxPyPzGlo(pVecNegVtx);
                  df4.getTrack(2).getPxPyPzGlo(pVecThirdVtx);
                  df4.getTrack(3).getPxPyPzGlo(pVecFourthVtx);
                  auto isBhadSel = helper.isSelectedBhadronToJPsi<4>(std::array{pVecPosVtx, pVecNegVtx, pVecThirdVtx, pVecFourthVtx}, std::array{track, trackFourth}, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexBto4tracks[0], secondaryVertexBto4tracks[1], secondaryVertexBto4tracks[2]}, activateQA, hMassVsPtB);
                  if (TESTBIT(isBhadSel, kB0ToJPsi)) {
                    keepEvent[kBtoJPsiKstar] = true;
                  }
                  if (TESTBIT(isBhadSel, kBsToJPsi)) {
                    keepEvent[kBtoJPsiPhi] = true;
                  }
                  if (TESTBIT(isBhadSel, kLbToJPsi)) {
                    keepEvent[kBtoJPsiPrKa] = true;
                  }
                }
              }
            }
          }
        } // end loop over tracks

        // 2-prong with Gamma (conversion photon)
        if (!keepEvent[kPhotonCharm2P] && isD0SignalTagged && (TESTBIT(selD0InMass, 0) || TESTBIT(selD0InMass, 1))) {
          auto photonsThisCollision = photons.sliceBy(photonsPerCollision, thisCollId);
          for (const auto& photon : photonsThisCollision) {
            auto posTrack = photon.posTrack_as<aod::V0Legs>();
            auto negTrack = photon.negTrack_as<aod::V0Legs>();
            if (!helper.isSelectedPhoton(photon, std::array{posTrack, negTrack}, activateQA, hV0Selected, hArmPod)) {
              continue;
            }
            std::array<float, 2> dcaInfo;
            std::array<float, 3> pVecPhoton = {photon.px(), photon.py(), photon.pz()};
            std::array<float, 3> posVecPhoton = {photon.vx(), photon.vy(), photon.vz()};
            auto trackParPhoton = o2::track::TrackPar(posVecPhoton, pVecPhoton, 0, true);
            trackParPhoton.setPID(o2::track::PID::Photon);
            trackParPhoton.setAbsCharge(0);
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPhoton, 2.f, matCorr, &dcaInfo);
            getPxPyPz(trackParPhoton, pVecPhoton);
            float massDStarCand{-1.}, massDStarBarCand{-999.};
            float massDiffDstar{-1.}, massDiffDstarBar{-999.};
            auto pVecReso2Prong = RecoDecay::pVec(pVec2Prong, pVecPhoton);
            auto ptCand = RecoDecay::pt(pVecReso2Prong);
            if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 1u)) {
              if (TESTBIT(selD0InMass, 0)) {
                massDStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecPhoton}, std::array{massPi, massKa, massGamma});
                massDiffDstar = massDStarCand - massD0Cand;
              }
              if (TESTBIT(selD0InMass, 1)) {
                massDStarBarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecPhoton}, std::array{massKa, massPi, massGamma});
                massDiffDstarBar = massDStarBarCand - massD0BarCand;
              }
              bool isGoodDstar = (cutsPtDeltaMassCharmReso->get(0u, 1u) < massDiffDstar && massDiffDstar < cutsPtDeltaMassCharmReso->get(1u, 1u));
              bool isGoodDstarBar = (cutsPtDeltaMassCharmReso->get(0u, 1u) < massDiffDstarBar && massDiffDstarBar < cutsPtDeltaMassCharmReso->get(1u, 1u));

              if (isGoodDstar || isGoodDstarBar) {
                if (activateQA) {
                  if (isGoodDstar) {
                    hMassVsPtC[kNCharmParticles + 1]->Fill(ptCand, massDiffDstar);
                  }
                  if (isGoodDstarBar) {
                    hMassVsPtC[kNCharmParticles + 1]->Fill(ptCand, massDiffDstarBar);
                  }
                }
                keepEvent[kPhotonCharm2P] = true;
                break; // we stop after the first D0-photon pair found
              }
            }
          }
        }

        // 2-prong with K0S or Lambda
        if (!keepEvent[kV0Charm2P] && isD0SignalTagged && (TESTBIT(selD0InMass, 0) || TESTBIT(selD0InMass, 1))) {
          auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, thisCollId);
          for (const auto& v0 : v0sThisCollision) {
            V0Cand v0Cand;
            if (!helper.buildV0(v0, tracksIU, collision, dfStrangeness, std::vector{cand2Prong.prong0Id(), cand2Prong.prong1Id()}, v0Cand)) {
              continue;
            }
            auto selV0 = helper.isSelectedV0(v0Cand, activateQA, hV0Selected, hArmPod);
            if (!selV0) {
              continue;
            }

            if (!keepEvent[kV0Charm2P] && TESTBIT(selV0, kK0S)) {

              // we first look for a D*+
              for (const auto& trackBachelorId : trackIdsThisCollision) { // start loop over tracks
                auto trackBachelor = tracks.rawIteratorAt(trackBachelorId.trackId());
                if (trackBachelor.globalIndex() == trackPos.globalIndex() || trackBachelor.globalIndex() == trackNeg.globalIndex() || trackBachelor.globalIndex() == v0.posTrackId() || trackBachelor.globalIndex() == v0.negTrackId()) {
                  continue;
                }

                auto trackParBachelor = getTrackPar(trackBachelor);
                std::array<float, 2> dcaBachelor{trackBachelor.dcaXY(), trackBachelor.dcaZ()};
                std::array<float, 3> pVecBachelor = trackBachelor.pVector();
                if (trackBachelor.collisionId() != thisCollId) {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelor, 2.f, noMatCorr, &dcaBachelor);
                  getPxPyPz(trackParBachelor, pVecBachelor);
                }

                auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty<kV0Charm2P>(trackBachelor, trackParBachelor, dcaBachelor);
                if (TESTBIT(isTrackSelected, kSoftPion) && ((TESTBIT(selD0InMass, 0) && trackBachelor.sign() > 0) || (TESTBIT(selD0InMass, 1) && trackBachelor.sign() < 0))) {
                  std::array<float, 2> massDausD0{massPi, massKa};
                  auto massD0dau = massD0Cand;
                  if (trackBachelor.sign() < 0) {
                    massDausD0[0] = massKa;
                    massDausD0[1] = massPi;
                    massD0dau = massD0BarCand;
                  }
                  auto pVecDStarCand = RecoDecay::pVec(pVec2Prong, pVecBachelor);
                  auto ptDStarCand = RecoDecay::pt(pVecDStarCand);
                  double massDStarCand{-999.}, massDiffDstar{-999.};
                  if (ptDStarCand > cutsPtDeltaMassCharmReso->get(2u, 0u)) {
                    massDStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor}, std::array{massDausD0[0], massDausD0[1], massPi});
                    massDiffDstar = massDStarCand - massD0dau;
                    if (cutsPtDeltaMassCharmReso->get(0u, 0u) <= massDiffDstar && massDiffDstar <= cutsPtDeltaMassCharmReso->get(1u, 0u)) {
                      if (activateQA) {
                        hMassVsPtC[kNCharmParticles]->Fill(ptDStarCand, massDiffDstar);
                      }
                      auto pVecReso2Prong = RecoDecay::pVec(pVecDStarCand, v0Cand.mom);
                      auto ptCand = RecoDecay::pt(pVecReso2Prong);
                      if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 3u)) {
                        auto massDStarK0S = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor, v0Cand.mom}, std::array{massDausD0[0], massDausD0[1], massPi, massK0S});
                        auto massDiffDsReso = massDStarK0S - massDStarCand;
                        if (cutsPtDeltaMassCharmReso->get(0u, 3u) < massDiffDsReso && massDiffDsReso < cutsPtDeltaMassCharmReso->get(1u, 3u)) {
                          if (activateQA) {
                            hMassVsPtC[kNCharmParticles + 3]->Fill(ptCand, massDiffDsReso);
                          }
                          keepEvent[kV0Charm2P] = true;
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
            if (!keepEvent[kV0Charm2P] && (TESTBIT(selV0, kLambda) || TESTBIT(selV0, kAntiLambda))) { // Xic(3055) and Xic(3080) --> since it occupies only a small bandwidth, we might want to keep also wrong sign pairs
              float massXicStarCand{-999.}, massXicStarBarCand{-999.};
              float massDiffXicStarCand{-999.}, massDiffXicStarBarCand{-999.};
              bool isRightSignXicStar{false}, isRightSignXicStarBar{false};
              auto pVecReso2Prong = RecoDecay::pVec(pVec2Prong, v0Cand.mom);
              auto ptCand = RecoDecay::pt(pVecReso2Prong);
              if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 5u)) {
                if (TESTBIT(selD0InMass, 0)) {
                  massXicStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, v0Cand.mom}, std::array{massPi, massKa, massLambda});
                  massDiffXicStarCand = massXicStarCand - massD0Cand;
                  isRightSignXicStar = TESTBIT(selV0, kLambda); // right sign if Lambda
                }
                if (TESTBIT(selD0InMass, 1)) {
                  massXicStarBarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, v0Cand.mom}, std::array{massKa, massPi, massLambda});
                  massDiffXicStarBarCand = massXicStarBarCand - massD0BarCand;
                  isRightSignXicStarBar = TESTBIT(selV0, kAntiLambda); // right sign if AntiLambda
                }
                bool isGoodXicStar = (cutsPtDeltaMassCharmReso->get(0u, 5u) < massDiffXicStarCand && massDiffXicStarCand < cutsPtDeltaMassCharmReso->get(1u, 5u));
                bool isGoodXicStarBar = (cutsPtDeltaMassCharmReso->get(0u, 5u) < massDiffXicStarBarCand && massDiffXicStarBarCand < cutsPtDeltaMassCharmReso->get(1u, 5u));

                if (activateQA) {
                  if (isGoodXicStar) {
                    if (isRightSignXicStar) {
                      hMassVsPtC[kNCharmParticles + 7]->Fill(ptCand, massDiffXicStarCand);
                    } else if (keepAlsoWrongDmesLambdaPairs) {
                      hMassVsPtC[kNCharmParticles + 8]->Fill(ptCand, massDiffXicStarBarCand);
                    }
                  }
                  if (isGoodXicStarBar) {
                    if (isRightSignXicStarBar) {
                      hMassVsPtC[kNCharmParticles + 7]->Fill(ptCand, massDiffXicStarCand);
                    } else if (keepAlsoWrongDmesLambdaPairs) {
                      hMassVsPtC[kNCharmParticles + 8]->Fill(ptCand, massDiffXicStarBarCand);
                    }
                  }
                }
                if ((isGoodXicStar && (isRightSignXicStar || keepAlsoWrongDmesLambdaPairs)) || (isGoodXicStarBar && (isRightSignXicStarBar || keepAlsoWrongDmesLambdaPairs))) {
                  keepEvent[kV0Charm2P] = true;
                  break;
                }
              }
            }
          }
        } // end V0 selection

        // 2-prong (D0 or D*) with proton for Lc resonances and ThetaC (3100)
        if (!keepEvent[kPrCharm2P] && isD0SignalTagged && (TESTBIT(selD0InMass, 0) || TESTBIT(selD0InMass, 1))) {
          for (const auto& trackProtonId : trackIdsThisCollision) { // start loop over tracks selecting only protons
            auto trackProton = tracks.rawIteratorAt(trackProtonId.trackId());
            auto trackParBachelorProton = getTrackPar(trackProton);
            if (trackProton.globalIndex() == trackPos.globalIndex() || trackProton.globalIndex() == trackNeg.globalIndex()) {
              continue;
            }
            std::array<float, 2> dcaInfoBachProton;
            if (trackProton.collisionId() != thisCollId) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelorProton, 2.f, noMatCorr, &dcaInfoBachProton);
            }
            std::array<float, 3> pVecProton = trackProton.pVector();
            bool isSelPIDProton = helper.isSelectedProton4CharmOrBeautyBaryons<false>(trackProton);
            if (isSelPIDProton) {
              if (!keepEvent[kPrCharm2P]) {
                // we first look for a D*+
                for (const auto& trackBachelorId : trackIdsThisCollision) { // start loop over tracks to find bachelor pion
                  if (!helper.isSelectedProtonFromLcResoOrThetaC<true>(trackProton)) {
                    continue;
                  } // stop here if proton below pT threshold for thetaC to avoid computational losses
                  auto trackBachelor = tracks.rawIteratorAt(trackBachelorId.trackId());
                  if (trackBachelor.globalIndex() == trackPos.globalIndex() || trackBachelor.globalIndex() == trackNeg.globalIndex() || trackBachelor.globalIndex() == trackProton.globalIndex()) {
                    continue;
                  }
                  auto trackParBachelor = getTrackPar(trackBachelor);
                  std::array<float, 2> dcaBachelor{trackBachelor.dcaXY(), trackBachelor.dcaZ()};
                  std::array<float, 3> pVecBachelor = trackBachelor.pVector();
                  if (trackBachelor.collisionId() != thisCollId) {
                    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelor, 2.f, noMatCorr, &dcaBachelor);
                    getPxPyPz(trackParBachelor, pVecBachelor);
                  }
                  auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty<kPrCharm2P>(trackBachelor, trackParBachelor, dcaBachelor);
                  if (TESTBIT(isTrackSelected, kSoftPion) && ((TESTBIT(selD0InMass, 0) && trackBachelor.sign() > 0) || (TESTBIT(selD0InMass, 1) && trackBachelor.sign() < 0))) {
                    if (pt2Prong < cutsPtDeltaMassCharmReso->get(3u, 12u)) {
                      continue;
                    }
                    std::array<float, 2> massDausD0{massPi, massKa};
                    auto massD0dau = massD0Cand;
                    if (trackBachelor.sign() < 0) {
                      massDausD0[0] = massKa;
                      massDausD0[1] = massPi;
                      massD0dau = massD0BarCand;
                    }
                    auto pVecDStarCand = RecoDecay::pVec(pVec2Prong, pVecBachelor);
                    auto ptDStarCand = RecoDecay::pt(pVecDStarCand);
                    double massDStarCand{-999.}, massDiffDstar{-999.};
                    if (ptDStarCand > cutsPtDeltaMassCharmReso->get(2u, 0u)) {
                      massDStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor}, std::array{massDausD0[0], massDausD0[1], massPi});
                      massDiffDstar = massDStarCand - massD0dau;
                      if (cutsPtDeltaMassCharmReso->get(0u, 0u) <= massDiffDstar && massDiffDstar <= cutsPtDeltaMassCharmReso->get(1u, 0u)) {
                        if (activateQA) { // probably this is not needed, since already performed for the Xic (duplicate)
                          hMassVsPtC[kNCharmParticles]->Fill(ptDStarCand, massDiffDstar);
                        }
                        auto pVecReso2Prong = RecoDecay::pVec(pVecDStarCand, pVecProton);
                        auto ptCand = RecoDecay::pt(pVecReso2Prong);
                        if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 12u)) {
                          // build D*0p candidate with the possibility of storing also the other sign hyp.
                          float massThetacCand{-999.}, massThetacBarCand{-999.};
                          float massDiffThetacCand{-999.}, massDiffThetacBarCand{-999.};
                          bool isRightSignThetaC{false}, isRightSignThetaCBar{false};
                          if (TESTBIT(selD0InMass, 1)) { // Correct hyp: ThetaC -> pD*- -> D0bar\pi-   (Equivalent to trackBachelor.sign() < 0)
                            massThetacCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor, pVecProton}, std::array{massDausD0[0], massDausD0[1], massPi, massProton});
                            massDiffThetacCand = massThetacCand - massDStarCand;
                            isRightSignThetaC = trackProton.sign() > 0; // right sign if proton
                          }
                          if (TESTBIT(selD0InMass, 0)) { // Correct hyp: ThetaCbar -> pD*+ -> pD0\pi+  (Equivalent to trackBachelor.sign() > 0)
                            massThetacBarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecBachelor, pVecProton}, std::array{massDausD0[0], massDausD0[1], massPi, massProton});
                            massDiffThetacBarCand = massThetacBarCand - massDStarCand;
                            isRightSignThetaCBar = trackProton.sign() < 0; // right sign if antiproton
                          }
                          bool isGoodThetac = (cutsPtDeltaMassCharmReso->get(0u, 12u) < massDiffThetacCand && massDiffThetacCand < cutsPtDeltaMassCharmReso->get(1u, 12u));
                          bool isGoodThetacBar = (cutsPtDeltaMassCharmReso->get(0u, 12u) < massDiffThetacBarCand && massDiffThetacBarCand < cutsPtDeltaMassCharmReso->get(1u, 12u));

                          if (activateQA) {
                            if (isGoodThetac) {
                              if (isRightSignThetaC) {
                                hMassVsPtC[kNCharmParticles + 21]->Fill(ptCand, massDiffThetacCand);
                              } else if (keepAlsoWrongDmesProtonPairs) {
                                hMassVsPtC[kNCharmParticles + 22]->Fill(ptCand, massDiffThetacBarCand);
                              }
                            }
                            if (isGoodThetacBar) {
                              if (isRightSignThetaCBar) {
                                hMassVsPtC[kNCharmParticles + 21]->Fill(ptCand, massDiffThetacCand);
                              } else if (keepAlsoWrongDmesProtonPairs) {
                                hMassVsPtC[kNCharmParticles + 22]->Fill(ptCand, massDiffThetacBarCand);
                              }
                            }
                          }
                          if ((isGoodThetac && (isRightSignThetaC || keepAlsoWrongDstarMesProtonPairs)) || (isGoodThetacBar && (isRightSignThetaCBar || keepAlsoWrongDstarMesProtonPairs))) {
                            keepEvent[kPrCharm2P] = true;
                            break;
                          }
                        }
                      }
                    }
                  }
                } // end bachelor pion for D*p pairs
                // build D0p candidate with the possibility of storing also the other sign hyp.
                if (pt2Prong < cutsPtDeltaMassCharmReso->get(3u, 11u)) {
                  continue;
                }
                if (!helper.isSelectedProtonFromLcResoOrThetaC(trackProton)) {
                  continue;
                }
                float massLcStarCand{-999.}, massLcStarBarCand{-999.};
                float massDiffLcStarCand{-999.}, massDiffLcStarBarCand{-999.};
                bool isRightSignLcStar{false}, isRightSignLcStarBar{false};
                auto pVecReso2Prong = RecoDecay::pVec(pVec2Prong, pVecProton);
                auto ptCand = RecoDecay::pt(pVecReso2Prong);
                if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 11u)) {
                  if (TESTBIT(selD0InMass, 0)) {
                    massLcStarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecProton}, std::array{massPi, massKa, massProton});
                    massDiffLcStarCand = massLcStarCand - massD0Cand;
                    isRightSignLcStar = trackProton.sign() > 0; // right sign if proton
                  }
                  if (TESTBIT(selD0InMass, 1)) {
                    massLcStarBarCand = RecoDecay::m(std::array{pVecPos, pVecNeg, pVecProton}, std::array{massKa, massPi, massProton});
                    massDiffLcStarBarCand = massLcStarBarCand - massD0BarCand;
                    isRightSignLcStarBar = trackProton.sign() < 0; // right sign if antiproton
                  }
                  bool isGoodLcStar = (cutsPtDeltaMassCharmReso->get(0u, 11u) < massDiffLcStarCand && massDiffLcStarCand < cutsPtDeltaMassCharmReso->get(1u, 11u));
                  bool isGoodLcStarBar = (cutsPtDeltaMassCharmReso->get(0u, 11u) < massDiffLcStarBarCand && massDiffLcStarBarCand < cutsPtDeltaMassCharmReso->get(1u, 11u));

                  if (activateQA) {
                    if (isGoodLcStar) {
                      if (isRightSignLcStar) {
                        hMassVsPtC[kNCharmParticles + 19]->Fill(ptCand, massDiffLcStarCand);
                      } else if (keepAlsoWrongDmesProtonPairs) {
                        hMassVsPtC[kNCharmParticles + 20]->Fill(ptCand, massDiffLcStarBarCand);
                      }
                    }
                    if (isGoodLcStarBar) {
                      if (isRightSignLcStarBar) {
                        hMassVsPtC[kNCharmParticles + 19]->Fill(ptCand, massDiffLcStarCand);
                      } else if (keepAlsoWrongDmesProtonPairs) {
                        hMassVsPtC[kNCharmParticles + 20]->Fill(ptCand, massDiffLcStarBarCand);
                      }
                    }
                  }
                  if ((isGoodLcStar && (isRightSignLcStar || keepAlsoWrongDmesProtonPairs)) || (isGoodLcStarBar && (isRightSignLcStarBar || keepAlsoWrongDmesProtonPairs))) {
                    keepEvent[kPrCharm2P] = true;
                    break;
                  }
                }
              }
            } // end proton loop
          }
        } // end Lc resonances via D0-proton decays

      } // end loop over 2-prong candidates

      std::vector<std::vector<int64_t>> indicesDau3Prong{}, indicesDau3ProngPrompt{};
      auto cand3ProngsThisColl = cand3Prongs.sliceBy(hf3ProngPerCollision, thisCollId);
      for (const auto& cand3Prong : cand3ProngsThisColl) { // start loop over 3 prongs
        std::array<int8_t, kNCharmParticles - 1> is3Prong = {
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DsToKKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::LcToPKPi),
          TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::XicToPKPi)};
        if (!std::accumulate(is3Prong.begin(), is3Prong.end(), 0)) { // check if it's a D+, Ds+, Lc+ or Xic+
          continue;
        }

        auto trackFirst = tracks.rawIteratorAt(cand3Prong.prong0Id());
        auto trackSecond = tracks.rawIteratorAt(cand3Prong.prong1Id());
        auto trackThird = tracks.rawIteratorAt(cand3Prong.prong2Id());

        auto trackParFirst = getTrackParCov(trackFirst);
        auto trackParSecond = getTrackParCov(trackSecond);
        auto trackParThird = getTrackParCov(trackThird);
        std::array<float, 2> dcaFirst{trackFirst.dcaXY(), trackFirst.dcaZ()};
        std::array<float, 2> dcaSecond{trackSecond.dcaXY(), trackSecond.dcaZ()};
        std::array<float, 2> dcaThird{trackThird.dcaXY(), trackThird.dcaZ()};
        std::array<float, 3> pVecFirst = trackFirst.pVector();
        std::array<float, 3> pVecSecond = trackSecond.pVector();
        std::array<float, 3> pVecThird = trackThird.pVector();
        if (trackFirst.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFirst, 2.f, noMatCorr, &dcaFirst);
          getPxPyPz(trackParFirst, pVecFirst);
        }
        if (trackSecond.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSecond, 2.f, noMatCorr, &dcaSecond);
          getPxPyPz(trackParSecond, pVecSecond);
        }
        if (trackThird.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
          getPxPyPz(trackParThird, pVecThird);
        }

        if (is3Prong[0]) { // D+ preselections
          is3Prong[0] = helper.isDplusPreselected(trackSecond);
        }
        if (is3Prong[1]) { // Ds preselections
          is3Prong[1] = helper.isDsPreselected(pVecFirst, pVecThird, pVecSecond, trackSecond);
        }
        if (is3Prong[2] || is3Prong[3]) { // charm baryon preselections
          auto presel = helper.isCharmBaryonPreselected(trackFirst, trackThird, trackSecond);
          if (is3Prong[2]) {
            is3Prong[2] = presel;
          }
          if (is3Prong[3]) {
            is3Prong[3] = presel;
          }
        }

        std::array<int8_t, kNCharmParticles - 1> isSignalTagged = is3Prong;
        std::array<int8_t, kNCharmParticles - 1> isCharmTagged = is3Prong;
        std::array<int8_t, kNCharmParticles - 1> isBeautyTagged = is3Prong;

        std::array<std::vector<float>, kNCharmParticles - 1> scores{};
        scores[0].insert(scores[0].end(), cand3Prong.mlProbSkimDplusToPiKPi().begin(), cand3Prong.mlProbSkimDplusToPiKPi().end());
        scores[1].insert(scores[1].end(), cand3Prong.mlProbSkimDsToKKPi().begin(), cand3Prong.mlProbSkimDsToKKPi().end());
        scores[2].insert(scores[2].end(), cand3Prong.mlProbSkimLcToPKPi().begin(), cand3Prong.mlProbSkimLcToPKPi().end());
        scores[3].insert(scores[3].end(), cand3Prong.mlProbSkimXicToPKPi().begin(), cand3Prong.mlProbSkimXicToPKPi().end());

        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles - 1; ++iCharmPart) {
          if (!is3Prong[iCharmPart]) { // we immediately skip if it was not selected for a given 3-prong species
            continue;
          }

          if (scores[iCharmPart].size() != 3) {
            scores[iCharmPart].resize(3);
            scores[iCharmPart][0] = 2.;
            scores[iCharmPart][1] = -1.;
            scores[iCharmPart][2] = -1.;
          }
          auto tagBDT = helper.isBDTSelected(scores[iCharmPart], thresholdBDTScores[iCharmPart + 1]);

          isCharmTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
          isBeautyTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);
          isSignalTagged[iCharmPart] = acceptBdtBkgOnly ? TESTBIT(tagBDT, RecoDecay::OriginType::None) : (isCharmTagged[iCharmPart] || isBeautyTagged[iCharmPart]);

          if (activateQA > 1) {
            hBDTScoreBkg[iCharmPart + 1]->Fill(scores[iCharmPart][0]);
            hBDTScorePrompt[iCharmPart + 1]->Fill(scores[iCharmPart][1]);
            hBDTScoreNonPrompt[iCharmPart + 1]->Fill(scores[iCharmPart][2]);
          }
        }

        if (!std::accumulate(isSignalTagged.begin(), isSignalTagged.end(), 0)) {
          continue;
        }

        keepEvent[kSingleCharm3P] = true;
        if (std::accumulate(isBeautyTagged.begin(), isBeautyTagged.end(), 0)) {
          keepEvent[kSingleNonPromptCharm3P] = true;
        }

        if (!keepOnlyDplusForDouble3Prongs) {
          indicesDau3Prong.push_back(std::vector<int64_t>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
          if (std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0)) {
            indicesDau3ProngPrompt.push_back(std::vector<int64_t>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
          }
        } else {
          if (isSignalTagged[kDplus - 1]) {
            indicesDau3Prong.push_back(std::vector<int64_t>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
            if (isCharmTagged[kDplus - 1]) {
              indicesDau3ProngPrompt.push_back(std::vector<int64_t>{trackFirst.globalIndex(), trackSecond.globalIndex(), trackThird.globalIndex()});
            }
          }
        } // end multiple 3-prong selection

        auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
        auto pt3Prong = RecoDecay::pt(pVec3Prong);
        float sign3Prong = -1 * trackFirst.sign() * trackSecond.sign() * trackThird.sign();

        std::array<int8_t, kNCharmParticles - 1> is3ProngInMass{0};
        if (is3Prong[0]) {
          is3ProngInMass[0] = helper.isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, activateQA, hMassVsPtC[kDplus]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, o2::constants::physics::Pdg::kDPlus, pt3Prong, scores[0][0], scores[0][1], scores[0][2]);
          }
        }
        if (is3Prong[1]) {
          is3ProngInMass[1] = helper.isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[1], activateQA, hMassVsPtC[kDs]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, o2::constants::physics::Pdg::kDS, pt3Prong, scores[1][0], scores[1][1], scores[1][2]);
          }
        }
        if (is3Prong[2]) {
          is3ProngInMass[2] = helper.isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[2], activateQA, hMassVsPtC[kLc]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, o2::constants::physics::Pdg::kLambdaCPlus, pt3Prong, scores[2][0], scores[2][1], scores[2][2]);
          }
        }
        if (is3Prong[3]) {
          is3ProngInMass[3] = helper.isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[3], activateQA, hMassVsPtC[kXic]);
          if (applyOptimisation) {
            optimisationTreeCharm(thisCollId, o2::constants::physics::Pdg::kXiCPlus, pt3Prong, scores[3][0], scores[3][1], scores[3][2]);
          }
        }

        if (helper.isSelectedHighPt3Prong(pt3Prong)) {
          keepEvent[kHighPt3P] = true;
          if (activateQA) {
            for (auto iCharmPart{1}; iCharmPart < kNCharmParticles; ++iCharmPart) {
              if (is3Prong[iCharmPart - 1] && (isSignalTagged[iCharmPart - 1])) {
                hCharmHighPt[iCharmPart]->Fill(pt3Prong);
              }
            }
          }
        } // end high-pT selection

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        auto tracksWithItsPid = soa::Attach<BigTracksPID, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

        for (const auto& trackId : trackIdsThisCollision) { // start loop over track indices as associated to this collision in HF code
          auto track = tracksWithItsPid.rawIteratorAt(trackId.trackId());
          if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
            continue;
          }

          auto trackParFourth = getTrackParCov(track);
          std::array<float, 2> dcaFourth{track.dcaXY(), track.dcaZ()};
          std::array<float, 3> pVecFourth = track.pVector();
          if (track.collisionId() != thisCollId) {
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFourth, 2.f, noMatCorr, &dcaFourth);
            getPxPyPz(trackParFourth, pVecFourth);
          }

          int charmParticleID[kNBeautyParticles - 3] = {o2::constants::physics::Pdg::kDPlus, o2::constants::physics::Pdg::kDS, o2::constants::physics::Pdg::kLambdaCPlus, o2::constants::physics::Pdg::kXiCPlus};

          float massCharmHypos[kNBeautyParticles - 3] = {massDPlus, massDs, massLc, massXic};
          auto isTrackSelected = helper.isSelectedTrackForSoftPionOrBeauty<kBeauty4P>(track, trackParFourth, dcaFourth);
          if (track.sign() * sign3Prong < 0 && TESTBIT(isTrackSelected, kForBeauty)) {
            for (int iHypo{0}; iHypo < kNBeautyParticles - 3 && !keepEvent[kBeauty4P]; ++iHypo) {
              if (isBeautyTagged[iHypo] && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1))) {
                auto massCandB = RecoDecay::m(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
                auto pVecBeauty4Prong = RecoDecay::pVec(pVec3Prong, pVecFourth);
                auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                if (helper.isSelectedBhadronInMassRange(ptCandBeauty4Prong, massCandB, iHypo + 3)) { // + 3 to account for B+ and B0->D*+ and Bc
                  if (activateQA) {
                    registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::Skimmed, iHypo + 3);
                  }
                  if (!activateSecVtxForB) {
                    keepEvent[kBeauty4P] = true;
                    if (applyOptimisation) {
                      optimisationTreeBeauty(thisCollId, charmParticleID[iHypo], pt3Prong, scores[iHypo][0], scores[iHypo][1], scores[iHypo][2], dcaFourth[0]);
                    }
                    if (activateQA) {
                      hMassVsPtB[iHypo + 3]->Fill(ptCandBeauty4Prong, massCandB);
                    }
                  } else {
                    int nVtxD{0};
                    try {
                      nVtxD = df3.process(trackParFirst, trackParSecond, trackParThird);
                    } catch (...) {
                      LOG(error) << "Exception caught in DCA fitter process call for charm 3-prong!";
                      nVtxD = 0;
                    }
                    if (nVtxD != 0) {
                      std::array<float, 3> pVecFirstVtx{}, pVecSecondVtx{}, pVecThirdVtx{};
                      df3.getTrack(0).getPxPyPzGlo(pVecFirstVtx);
                      df3.getTrack(1).getPxPyPzGlo(pVecSecondVtx);
                      df3.getTrack(1).getPxPyPzGlo(pVecThirdVtx);
                      auto trackParD = df3.createParentTrackParCov();
                      trackParD.setAbsCharge(sign3Prong); // to be sure
                      auto pVec3ProngVtx = RecoDecay::pVec(pVecFirstVtx, pVecSecondVtx, pVecThirdVtx);
                      int nVtxB{0};
                      try {
                        nVtxB = dfB.process(trackParD, trackParFourth);
                      } catch (...) {
                        LOG(error) << "Exception caught in DCA fitter process call for B 4-prong!";
                        nVtxB = 0;
                      }
                      if (nVtxB != 0) {
                        if (activateQA) {
                          registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::BeautyVertex, iHypo + 3);
                        }
                        const auto& secondaryVertexB = dfB.getPCACandidate();
                        std::array<float, 3> pVecFourtVtx{};
                        dfB.getTrack(0).getPxPyPzGlo(pVec3ProngVtx);
                        dfB.getTrack(1).getPxPyPzGlo(pVecFourtVtx);
                        std::array<float, 2> dca3Prong; //{trackParD.dcaXY(), trackParD.dcaZ()};
                        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParD, 2.f, noMatCorr, &dca3Prong);
                        bool isBhad = helper.isSelectedBhadron(pVec3ProngVtx, pVecFourtVtx, dca3Prong, dcaFourth, std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexB[0], secondaryVertexB[1], secondaryVertexB[2]}, iHypo + 3);
                        if (isBhad) {
                          keepEvent[kBeauty4P] = true;
                          // fill optimisation tree
                          if (applyOptimisation) {
                            optimisationTreeBeauty(thisCollId, charmParticleID[iHypo], pt3Prong, scores[iHypo][0], scores[iHypo][1], scores[iHypo][2], dcaFourth[0]);
                          }
                          if (activateQA) {
                            registry.fill(HIST("fHfVtxStages"), 1 + HfVtxStage::CharmHadPiSelected, iHypo + 3);
                            hCpaVsPtB[iHypo + 3]->Fill(ptCandBeauty4Prong, RecoDecay::cpa(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexB[0], secondaryVertexB[1], secondaryVertexB[2]}, RecoDecay::pVec(pVec3ProngVtx, pVecFourtVtx)));
                            hDecayLengthVsPtB[iHypo + 3]->Fill(ptCandBeauty4Prong, RecoDecay::distance(std::array<double, 3>{static_cast<double>(collision.posX()), static_cast<double>(collision.posY()), static_cast<double>(collision.posZ())}, std::array{secondaryVertexB[0], secondaryVertexB[1], secondaryVertexB[2]}));
                            hImpactParamProductVsPtB[iHypo + 3]->Fill(ptCandBeauty4Prong, dca3Prong[0] * dcaFourth[0]);
                            hMassVsPtB[iHypo + 3]->Fill(ptCandBeauty4Prong, massCandB);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } // end beauty selection

          // 3-prong femto
          bool isProton = helper.isSelectedTrack4Femto(track, trackParFourth, activateQA, hPrDePID[0], hPrDePID[1], kProtonForFemto);
          bool isDeuteron = helper.isSelectedTrack4Femto(track, trackParFourth, activateQA, hPrDePID[2], hPrDePID[3], kDeuteronForFemto);

          if (isProton && track.collisionId() == thisCollId) {
            for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
              if (isCharmTagged[iHypo] && enableFemtoChannels->get(0u, iHypo + 1)) {
                float relativeMomentum = helper.computeRelativeMomentum(pVecFourth, pVec3Prong, massCharmHypos[iHypo]);
                if (applyOptimisation) {
                  optimisationTreeFemto(thisCollId, charmParticleID[iHypo], pt3Prong, scores[iHypo][0], scores[iHypo][1], scores[iHypo][2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr(), track.tpcNSigmaDe(), track.tofNSigmaDe());
                }
                if (relativeMomentum < femtoMaxRelativeMomentum) {
                  keepEvent[kFemto3P] = true;
                  if (activateQA) {
                    hCharmProtonKstarDistr[iHypo + 1]->Fill(relativeMomentum);
                  }
                }
              }
            }
          }
          if (isDeuteron && track.collisionId() == thisCollId) {
            for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
              if (isCharmTagged[iHypo] && enableFemtoChannels->get(1u, iHypo + 1)) {
                float relativeMomentum = helper.computeRelativeMomentum(pVecFourth, pVec3Prong, massCharmHypos[iHypo]);
                if (applyOptimisation) {
                  optimisationTreeFemto(thisCollId, charmParticleID[iHypo], pt3Prong, scores[iHypo][0], scores[iHypo][1], scores[iHypo][2], relativeMomentum, track.tpcNSigmaPr(), track.tofNSigmaPr(), track.tpcNSigmaDe(), track.tofNSigmaDe());
                }
                if (relativeMomentum < femtoMaxRelativeMomentum) {
                  keepEvent[kFemto3P] = true;
                  if (activateQA) {
                    hCharmDeuteronKstarDistr[iHypo + 1]->Fill(relativeMomentum);
                  }
                }
              }
            }
          } // end femto selection

          // SigmaC++ K- trigger
          if (!keepEvent[kSigmaCPPK] && is3Prong[2] > 0 && is3ProngInMass[2] > 0 && isSignalTagged[2] > 0 && helper.isSelectedKaonFromXicResoToSigmaC<true>(track)) {
            // we need a candidate Lc->pKpi and a candidate soft kaon

            // look for SigmaC++ candidates
            for (const auto& trackSoftPiId : trackIdsThisCollision) { // start loop over tracks (soft pi)

              // soft pion candidates
              auto trackSoftPi = tracks.rawIteratorAt(trackSoftPiId.trackId());
              auto globalIndexSoftPi = trackSoftPi.globalIndex();

              // exclude tracks already used to build the 3-prong candidate
              if (globalIndexSoftPi == trackFirst.globalIndex() || globalIndexSoftPi == trackSecond.globalIndex() || globalIndexSoftPi == trackThird.globalIndex()) {
                // do not consider as candidate soft pion a track already used to build the current 3-prong candidate
                continue;
              }

              // exclude already the current track if it corresponds to the K- candidate
              if (globalIndexSoftPi == track.globalIndex()) {
                continue;
              }

              // check the candidate SigmaC++ charge
              std::array<int, 4> chargesSc = {trackFirst.sign(), trackSecond.sign(), trackThird.sign(), trackSoftPi.sign()};
              int chargeSc = std::accumulate(chargesSc.begin(), chargesSc.end(), 0); // SIGNED electric charge of SigmaC candidate
              if (std::abs(chargeSc) != 2) {
                continue;
              }

              // select soft pion candidates
              auto trackParSoftPi = getTrackPar(trackSoftPi);
              std::array<float, 2> dcaSoftPi{trackSoftPi.dcaXY(), trackSoftPi.dcaZ()};
              std::array<float, 3> pVecSoftPi = trackSoftPi.pVector();
              if (trackSoftPi.collisionId() != thisCollId) {
                // This is a track reassociated to this PV by the track-to-collision-associator
                // Let's propagate this track to it, and calculate dcaXY, dcaZ
                o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSoftPi, 2.f, noMatCorr, &dcaSoftPi);
                getPxPyPz(trackParSoftPi, pVecSoftPi);
              }
              int16_t isSoftPionSelected = helper.isSelectedTrackForSoftPionOrBeauty<kSigmaCPPK>(trackSoftPi, trackParSoftPi, dcaSoftPi);
              if (TESTBIT(isSoftPionSelected, kSoftPionForSigmaC) /*&& (TESTBIT(is3Prong[2], 0) || TESTBIT(is3Prong[2], 1))*/) {

                // check the mass of the SigmaC++ candidate
                auto pVecSigmaC = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird, pVecSoftPi);
                auto ptSigmaC = RecoDecay::pt(pVecSigmaC);
                int8_t whichSigmaC = helper.isSelectedSigmaCInDeltaMassRange<2>(pVecFirst, pVecThird, pVecSecond, pVecSoftPi, ptSigmaC, is3Prong[2], hMassVsPtC[kNCharmParticles + 9], activateQA);
                if (whichSigmaC > 0) {
                  /// let's build a candidate SigmaC++K- pair
                  /// and keep it only if:
                  ///   - it has the correct charge (±1)
                  ///   - it is in the correct mass range

                  // check the charge for SigmaC++K- candidates
                  if (std::abs(chargeSc + track.sign()) != 1) {
                    continue;
                  }

                  // check the invariant mass
                  float massSigmaCPKPi{-999.}, massSigmaCPiKP{-999.}, deltaMassXicResoPKPi{-999.}, deltaMassXicResoPiKP{-999.};
                  float ptSigmaCKaon = RecoDecay::pt(pVecSigmaC, pVecFourth);

                  if (ptSigmaCKaon > cutsPtDeltaMassCharmReso->get(2u, 10u)) {
                    if (TESTBIT(whichSigmaC, 0)) {
                      massSigmaCPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi}, std::array{massProton, massKa, massPi, massPi});
                      deltaMassXicResoPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi, pVecFourth}, std::array{massProton, massKa, massPi, massPi, massKa}) - massSigmaCPKPi;
                    }
                    if (TESTBIT(whichSigmaC, 1)) {
                      massSigmaCPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi}, std::array{massPi, massKa, massProton, massPi});
                      deltaMassXicResoPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi, pVecFourth}, std::array{massPi, massKa, massProton, massPi, massKa}) - massSigmaCPiKP;
                    }
                    bool isPKPiOk = (cutsPtDeltaMassCharmReso->get(0u, 10u) < deltaMassXicResoPKPi && deltaMassXicResoPKPi < cutsPtDeltaMassCharmReso->get(1u, 10u));
                    bool isPiKPOk = (cutsPtDeltaMassCharmReso->get(0u, 10u) < deltaMassXicResoPiKP && deltaMassXicResoPiKP < cutsPtDeltaMassCharmReso->get(1u, 10u));
                    if (isPKPiOk || isPiKPOk) {
                      /// This is a good SigmaC++K- event
                      keepEvent[kSigmaCPPK] = true;

                      /// QA plot
                      if (activateQA) {
                        if (isPKPiOk) {
                          if (TESTBIT(whichSigmaC, 2)) {
                            hMassVsPtC[kNCharmParticles + 11]->Fill(ptSigmaCKaon, deltaMassXicResoPKPi);
                          }
                          if (TESTBIT(whichSigmaC, 3)) {
                            hMassVsPtC[kNCharmParticles + 12]->Fill(ptSigmaCKaon, deltaMassXicResoPKPi);
                          }
                        }
                        if (isPiKPOk) {
                          if (TESTBIT(whichSigmaC, 2)) {
                            hMassVsPtC[kNCharmParticles + 11]->Fill(ptSigmaCKaon, deltaMassXicResoPiKP);
                          }
                          if (TESTBIT(whichSigmaC, 3)) {
                            hMassVsPtC[kNCharmParticles + 12]->Fill(ptSigmaCKaon, deltaMassXicResoPiKP);
                          }
                        }
                      }
                    }
                  }
                }
              } // end SigmaC++ candidate
            } // end loop over tracks (soft pi)
          } // end candidate Lc->pKpi
        } // end loop over tracks

        // Ds with photon
        bool isGoodDsToKKPi = (isSignalTagged[kDs - 1]) && TESTBIT(is3ProngInMass[kDs - 1], 0);
        bool isGoodDsToPiKK = (isSignalTagged[kDs - 1]) && TESTBIT(is3ProngInMass[kDs - 1], 1);
        if (!keepEvent[kPhotonCharm3P] && (isGoodDsToKKPi || isGoodDsToPiKK)) {
          auto massDsKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massKa, massKa, massPi});
          auto massDsPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massKa});
          auto photonsThisCollision = photons.sliceBy(photonsPerCollision, thisCollId);
          for (const auto& photon : photonsThisCollision) {
            auto posTrack = photon.posTrack_as<aod::V0Legs>();
            auto negTrack = photon.negTrack_as<aod::V0Legs>();
            if (!helper.isSelectedPhoton(photon, std::array{posTrack, negTrack}, activateQA, hV0Selected, hArmPod)) {
              continue;
            }
            std::array<float, 2> dcaInfo;
            std::array<float, 3> pVecPhoton = {photon.px(), photon.py(), photon.pz()};
            std::array<float, 3> posVecPhoton = {photon.vx(), photon.vy(), photon.vz()};
            auto trackParPhoton = o2::track::TrackPar(posVecPhoton, pVecPhoton, 0, true);
            trackParPhoton.setAbsCharge(0);
            trackParPhoton.setPID(o2::track::PID::Photon);
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPhoton, 2.f, matCorr, &dcaInfo);
            getPxPyPz(trackParPhoton, pVecPhoton);
            float massDsStarToKKPiCand{-1.}, massDsStarToPiKKCand{999.};
            float massDiffDsStarToKKPi{-1.}, massDiffDsStarToPiKK{999.};
            if (isGoodDsToKKPi) {
              massDsStarToKKPiCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecPhoton}, std::array{massKa, massKa, massPi, massGamma});
              massDiffDsStarToKKPi = massDsStarToKKPiCand - massDsKKPi;
            }
            if (isGoodDsToPiKK) {
              massDsStarToPiKKCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecPhoton}, std::array{massPi, massKa, massKa, massGamma});
              massDiffDsStarToPiKK = massDsStarToPiKKCand - massDsPiKK;
            }

            auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, pVecPhoton);
            auto ptCand = RecoDecay::pt(pVecReso3Prong);
            if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 2u)) {
              bool isGoodDsStarToKKPi = (cutsPtDeltaMassCharmReso->get(0u, 2u) < massDiffDsStarToKKPi && massDiffDsStarToKKPi < cutsPtDeltaMassCharmReso->get(1u, 2u));
              bool isGoodDsStarToPiKK = (cutsPtDeltaMassCharmReso->get(0u, 2u) < massDiffDsStarToPiKK && massDiffDsStarToPiKK < cutsPtDeltaMassCharmReso->get(1u, 2u));
              if (isGoodDsStarToKKPi || isGoodDsStarToPiKK) {
                if (activateQA) {
                  if (isGoodDsStarToKKPi) {
                    hMassVsPtC[kNCharmParticles + 2]->Fill(ptCand, massDiffDsStarToKKPi);
                  }
                  if (isGoodDsStarToPiKK) {
                    hMassVsPtC[kNCharmParticles + 2]->Fill(ptCand, massDiffDsStarToPiKK);
                  }
                }
                keepEvent[kPhotonCharm3P] = true;
                break; // we stop after the first Ds + photon found
              }
            }
          }
        }

        // D+ with K0S or Lambda and SigmaC0 with K0S
        auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, thisCollId);
        bool isGoodDPlus = (isSignalTagged[kDplus - 1]) && is3ProngInMass[kDplus - 1];
        bool isGoodLcToPKPi = (isSignalTagged[kLc - 1]) && TESTBIT(is3ProngInMass[kLc - 1], 0);
        bool isGoodLcToPiKP = (isSignalTagged[kLc - 1]) && TESTBIT(is3ProngInMass[kLc - 1], 1);
        auto massDPlusCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massPi});

        if ((!keepEvent[kV0Charm3P] && isGoodDPlus) || (!keepEvent[kSigmaC0K0] && (isGoodLcToPKPi || isGoodLcToPiKP))) {
          for (const auto& v0 : v0sThisCollision) {
            V0Cand v0Cand;
            if (!helper.buildV0(v0, tracksIU, collision, dfStrangeness, std::vector{cand3Prong.prong0Id(), cand3Prong.prong1Id(), cand3Prong.prong2Id()}, v0Cand)) {
              continue;
            }
            auto selV0 = helper.isSelectedV0(v0Cand, activateQA, hV0Selected, hArmPod);
            if (!selV0) {
              continue;
            }

            // we pair D+ with V0
            if (!keepEvent[kV0Charm3P] && isGoodDPlus) {
              if (!keepEvent[kV0Charm3P] && TESTBIT(selV0, kK0S)) { // Ds2*
                auto massDsStarCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, v0Cand.mom}, std::array{massPi, massKa, massPi, massK0S});
                auto massDiffDsStar = massDsStarCand - massDPlusCand;
                auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, v0Cand.mom);
                auto ptCand = RecoDecay::pt(pVecReso3Prong);
                if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 4u)) {
                  if (cutsPtDeltaMassCharmReso->get(0u, 4u) < massDiffDsStar && massDiffDsStar < cutsPtDeltaMassCharmReso->get(1u, 4u)) {
                    if (activateQA) {
                      hMassVsPtC[kNCharmParticles + 4]->Fill(ptCand, massDiffDsStar);
                    }
                    keepEvent[kV0Charm3P] = true;
                  }
                }
              }
              if (!keepEvent[kV0Charm3P] && (TESTBIT(selV0, kLambda) || TESTBIT(selV0, kAntiLambda))) { // Xic(3055) and Xic(3080) --> since it occupies only a small bandwidth, we might want to keep also wrong sign pairs
                auto pVecReso3Prong = RecoDecay::pVec(pVec3Prong, v0Cand.mom);
                auto ptCand = RecoDecay::pt(pVecReso3Prong);
                if (ptCand > cutsPtDeltaMassCharmReso->get(2u, 5u)) {
                  auto massXicStarCand = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, v0Cand.mom}, std::array{massPi, massKa, massPi, massLambda});
                  auto massDiffXicStar = massXicStarCand - massDPlusCand;
                  bool isRightSign = ((TESTBIT(selV0, kLambda) && sign3Prong > 0) || (TESTBIT(selV0, kAntiLambda) && sign3Prong < 0));
                  if (cutsPtDeltaMassCharmReso->get(0u, 5u) < massDiffXicStar && massDiffXicStar < cutsPtDeltaMassCharmReso->get(1u, 5u)) {
                    if (activateQA) {
                      if (isRightSign) {
                        hMassVsPtC[kNCharmParticles + 5]->Fill(ptCand, massDiffXicStar);
                      } else if (keepAlsoWrongDmesLambdaPairs) {
                        hMassVsPtC[kNCharmParticles + 6]->Fill(ptCand, massDiffXicStar);
                      }
                    }
                    if (isRightSign || keepAlsoWrongDmesLambdaPairs) {
                      keepEvent[kV0Charm3P] = true;
                    }
                  }
                }
              }
            } // end D+ with V0

            // we pair SigmaC0 with V0
            if (!keepEvent[kSigmaC0K0] && (isGoodLcToPKPi || isGoodLcToPiKP) && TESTBIT(selV0, kK0S)) {
              // look for SigmaC0 candidates
              for (const auto& trackSoftPiId : trackIdsThisCollision) { // start loop over tracks (soft pi)

                // soft pion candidates
                auto trackSoftPi = tracks.rawIteratorAt(trackSoftPiId.trackId());
                auto globalIndexSoftPi = trackSoftPi.globalIndex();

                // exclude tracks already used to build the 3-prong candidate
                if (globalIndexSoftPi == trackFirst.globalIndex() || globalIndexSoftPi == trackSecond.globalIndex() || globalIndexSoftPi == trackThird.globalIndex() || globalIndexSoftPi == v0.posTrackId() || globalIndexSoftPi == v0.negTrackId()) {
                  // do not consider as candidate soft pion a track already used to build the current 3-prong candidate / V0 candidate
                  continue;
                }

                // check the candidate SigmaC0 charge
                std::array<int, 4> chargesSc = {trackFirst.sign(), trackSecond.sign(), trackThird.sign(), trackSoftPi.sign()};
                int chargeSc = std::accumulate(chargesSc.begin(), chargesSc.end(), 0); // SIGNED electric charge of SigmaC candidate
                if (chargeSc != 0) {
                  continue;
                }

                // select soft pion candidates
                auto trackParSoftPi = getTrackPar(trackSoftPi);
                std::array<float, 2> dcaSoftPi{trackSoftPi.dcaXY(), trackSoftPi.dcaZ()};
                std::array<float, 3> pVecSoftPi = trackSoftPi.pVector();
                if (trackSoftPi.collisionId() != thisCollId) {
                  // This is a track reassociated to this PV by the track-to-collision-associator
                  // Let's propagate this track to it, and calculate dcaXY, dcaZ
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSoftPi, 2.f, noMatCorr, &dcaSoftPi);
                  getPxPyPz(trackParSoftPi, pVecSoftPi);
                }
                int16_t isSoftPionSelected = helper.isSelectedTrackForSoftPionOrBeauty<kSigmaC0K0>(trackSoftPi, trackParSoftPi, dcaSoftPi);
                if (TESTBIT(isSoftPionSelected, kSoftPionForSigmaC) /*&& (TESTBIT(is3Prong[2], 0) || TESTBIT(is3Prong[2], 1))*/) {

                  // check the mass of the SigmaC0 candidate
                  auto pVecSigmaC = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird, pVecSoftPi);
                  auto ptSigmaC = RecoDecay::pt(pVecSigmaC);
                  int8_t whichSigmaC = helper.isSelectedSigmaCInDeltaMassRange<0>(pVecFirst, pVecThird, pVecSecond, pVecSoftPi, ptSigmaC, is3Prong[2], hMassVsPtC[kNCharmParticles + 10], activateQA);
                  if (whichSigmaC > 0) {
                    /// let's build a candidate SigmaC0K0s pair
                    /// and keep it only if it is in the correct mass range

                    float massSigmaCPKPi{-999.}, massSigmaCPiKP{-999.}, deltaMassXicResoPKPi{-999.}, deltaMassXicResoPiKP{-999.};
                    float ptSigmaCKaon = RecoDecay::pt(pVecSigmaC, v0Cand.mom);
                    if (ptSigmaCKaon > cutsPtDeltaMassCharmReso->get(2u, 10u)) {
                      if (TESTBIT(whichSigmaC, 0)) {
                        massSigmaCPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi}, std::array{massProton, massKa, massPi, massPi});
                        deltaMassXicResoPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi, v0Cand.mom}, std::array{massProton, massKa, massPi, massPi, massK0S}) - massSigmaCPKPi;
                      }
                      if (TESTBIT(whichSigmaC, 1)) {
                        massSigmaCPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi}, std::array{massPi, massKa, massProton, massPi});
                        deltaMassXicResoPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird, pVecSoftPi, v0Cand.mom}, std::array{massPi, massKa, massProton, massPi, massK0S}) - massSigmaCPiKP;
                      }

                      bool isPKPiOk = (cutsPtDeltaMassCharmReso->get(0u, 10u) < deltaMassXicResoPKPi && deltaMassXicResoPKPi < cutsPtDeltaMassCharmReso->get(1u, 10u));
                      bool isPiKPOk = (cutsPtDeltaMassCharmReso->get(0u, 10u) < deltaMassXicResoPiKP && deltaMassXicResoPiKP < cutsPtDeltaMassCharmReso->get(1u, 10u));
                      if (isPKPiOk || isPiKPOk) {
                        /// This is a good SigmaC0K0s event
                        keepEvent[kSigmaC0K0] = true;

                        /// QA plot
                        if (activateQA) {
                          if (isPKPiOk) {
                            if (TESTBIT(whichSigmaC, 2)) {
                              hMassVsPtC[kNCharmParticles + 13]->Fill(ptSigmaCKaon, deltaMassXicResoPKPi);
                            }
                            if (TESTBIT(whichSigmaC, 3)) {
                              hMassVsPtC[kNCharmParticles + 14]->Fill(ptSigmaCKaon, deltaMassXicResoPKPi);
                            }
                          }
                          if (isPiKPOk) {
                            if (TESTBIT(whichSigmaC, 2)) {
                              hMassVsPtC[kNCharmParticles + 13]->Fill(ptSigmaCKaon, deltaMassXicResoPiKP);
                            }
                            if (TESTBIT(whichSigmaC, 3)) {
                              hMassVsPtC[kNCharmParticles + 14]->Fill(ptSigmaCKaon, deltaMassXicResoPiKP);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } // end loop over tracks (soft pi)
            }
          }
        }
      } // end loop over 3-prong candidates

      if (!keepEvent[kCharmBarToXiBach] || !keepEvent[kCharmBarToXi2Bach]) {
        auto cascThisColl = cascades.sliceBy(cascPerCollision, thisCollId);
        for (const auto& casc : cascThisColl) {

          bool hasStrangeTrack{false};

          TracksIUPID::iterator cascTrack;
          int requireStrangenessTrackingAny = requireStrangenessTracking->get(0u, 0u) + requireStrangenessTracking->get(0u, 1u);
          if (requireStrangenessTrackingAny > 0) { // enabled for at least one of the two
            auto trackedCascIdThisColl = trackedCasc.sliceBy(trackedCascadesPerCollision, thisCollId);
            for (const auto& trackedCascId : trackedCascIdThisColl) {
              if (trackedCascId.cascadeId() == casc.globalIndex()) {
                hasStrangeTrack = true;
                cascTrack = trackedCascId.track_as<TracksIUPID>();
                break;
              }
            }
          }

          CascCand cascCand;
          if (!helper.buildCascade(casc, v0s, tracksIU, collision, dfStrangeness, {}, cascCand)) {
            continue;
          }

          if (!helper.isSelectedCascade(cascCand)) {
            continue;
          }

          if (activateQA) {
            hMassXi[0]->Fill(cascCand.sign, cascCand.mXi);
            if (hasStrangeTrack) {
              hMassXi[1]->Fill(cascCand.sign, cascCand.mXi);
            }
          }

          auto bachelorCascId = casc.bachelorId();
          auto v0 = v0s.rawIteratorAt(casc.v0Id());
          auto v0DauPosId = v0.posTrackId();
          auto v0DauNegId = v0.negTrackId();

          // propagate to PV
          std::array<float, 2> dcaInfo;
          o2::track::TrackParCov trackParCasc;
          o2::track::TrackParCov trackParCascTrack;
          if (requireStrangenessTrackingAny < 2) { // needed for at least one of the two
            trackParCasc = o2::track::TrackParCov(cascCand.vtx, cascCand.mom, cascCand.cov, cascCand.sign, true);
            trackParCasc.setPID(o2::track::PID::XiMinus);
            trackParCasc.setAbsCharge(1); // to be sure
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCasc, 2.f, matCorr, &dcaInfo);
          }
          if (requireStrangenessTrackingAny > 0 && hasStrangeTrack) { // needed for at least one of the two
            trackParCascTrack = getTrackParCov(cascTrack);
            trackParCascTrack.setPID(o2::track::PID::XiMinus);
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCascTrack, 2.f, matCorr, &dcaInfo);
          }

          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
          for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks (first bachelor)
            auto track = tracks.rawIteratorAt(trackId.trackId());

            // check if track is one of the Xi daughters
            if (track.globalIndex() == bachelorCascId || track.globalIndex() == v0DauPosId || track.globalIndex() == v0DauNegId) {
              continue;
            }

            auto trackParBachelor = getTrackParCov(track);
            std::array<float, 2> dcaInfoBach;
            if (track.collisionId() != thisCollId) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelor, 2.f, noMatCorr, &dcaInfoBach);
            }

            auto isSelBachelor = helper.isSelectedBachelorForCharmBaryon(track, dcaInfoBach);
            if (isSelBachelor == kRejected) {
              continue;
            }

            if (!keepEvent[kCharmBarToXiBach] && track.sign() * cascCand.sign < 0) { // XiPi and XiKa

              bool isSelXiBach{false};
              if (requireStrangenessTracking->get(0u, 0u) > 0) {
                if (hasStrangeTrack) {
                  isSelXiBach = helper.isSelectedXiBach(trackParCascTrack, trackParBachelor, isSelBachelor, collision, dfStrangeness, activateQA, hMassVsPtC[kNCharmParticles + 15], hMassVsPtC[kNCharmParticles + 16]);
                }
              } else {
                isSelXiBach = helper.isSelectedXiBach(trackParCasc, trackParBachelor, isSelBachelor, collision, dfStrangeness, activateQA, hMassVsPtC[kNCharmParticles + 15], hMassVsPtC[kNCharmParticles + 16]);
              }
              if (isSelXiBach) {
                keepEvent[kCharmBarToXiBach] = true;
              }
            }

            // only pions needed below
            if (!TESTBIT(isSelBachelor, kPionForCharmBaryon)) {
              continue;
            }

            if (!keepEvent[kCharmBarToXi2Bach]) {
              for (const auto& trackIdSecond : trackIdsThisCollision) { // start loop over tracks (second bachelor)
                auto trackSecond = tracks.rawIteratorAt(trackIdSecond.trackId());

                // check if track is one of the Xi daughters
                if (trackSecond.globalIndex() == track.globalIndex() || trackSecond.globalIndex() == bachelorCascId || trackSecond.globalIndex() == v0DauPosId || trackSecond.globalIndex() == v0DauNegId) {
                  continue;
                }

                if (track.sign() * trackSecond.sign() < 0 || track.sign() * cascCand.sign > 0) { // we want same sign pions, opposite to the xi
                  continue;
                }

                auto trackParBachelorSecond = getTrackParCov(trackSecond);
                std::array<float, 2> dcaInfoBachSecond;
                if (trackSecond.collisionId() != thisCollId) {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParBachelorSecond, 2.f, noMatCorr, &dcaInfoBachSecond);
                }

                auto isSelBachelorSecond = helper.isSelectedBachelorForCharmBaryon(trackSecond, dcaInfoBachSecond);
                if (!TESTBIT(isSelBachelorSecond, kPionForCharmBaryon)) {
                  continue;
                }
                if (!keepEvent[kCharmBarToXi2Bach]) { // XiPiPi

                  bool isSelXiBachBach{false};
                  if (requireStrangenessTracking->get(0u, 1u) > 0) {
                    if (hasStrangeTrack) {
                      isSelXiBachBach = helper.isSelectedXiBachBach<3>(trackParCascTrack, {trackParBachelor, trackParBachelorSecond}, collision, dfStrangeness3, activateQA, hMassVsPtC[kNCharmParticles + 17]);
                    }
                  } else { // vertex with only the two bachelors
                    isSelXiBachBach = helper.isSelectedXiBachBach<2>(trackParCasc, {trackParBachelor, trackParBachelorSecond}, collision, df2, activateQA, hMassVsPtC[kNCharmParticles + 17]);
                  }
                  if (isSelXiBachBach) {
                    keepEvent[kCharmBarToXi2Bach] = true;
                  }
                }
              }
            }
          }
        }
      }

      auto n2Prongs = helper.computeNumberOfCandidates(indicesDau2Prong);
      auto n2ProngsPrompt = helper.computeNumberOfCandidates(indicesDau2ProngPrompt);
      auto n3Prongs = helper.computeNumberOfCandidates(indicesDau3Prong);
      auto n3ProngsPrompt = helper.computeNumberOfCandidates(indicesDau3ProngPrompt);
      indicesDau2Prong.insert(indicesDau2Prong.end(), indicesDau3Prong.begin(), indicesDau3Prong.end());
      auto n23Prongs = helper.computeNumberOfCandidates(indicesDau2Prong);
      indicesDau2ProngPrompt.insert(indicesDau2ProngPrompt.end(), indicesDau3ProngPrompt.begin(), indicesDau3ProngPrompt.end());
      auto n23ProngsPrompt = helper.computeNumberOfCandidates(indicesDau2ProngPrompt);

      if (activateQA) {
        hN2ProngCharmCand->Fill(n2Prongs);
        hN3ProngCharmCand->Fill(n3Prongs);
      }

      if (n2Prongs > 1 && enableDoubleCharmChannels->get(0u, 0u)) {
        if (enableDoubleCharmChannels->get(1u, 0u)) {
          keepEvent[kDoubleCharm2P] = true;
        } else {
          if (n2ProngsPrompt > 1) {
            keepEvent[kDoubleCharm2P] = true;
          }
        }
      }
      if (n3Prongs > 1 && enableDoubleCharmChannels->get(0u, 1u)) {
        if (enableDoubleCharmChannels->get(1u, 1u)) {
          keepEvent[kDoubleCharm3P] = true;
        } else {
          if (n3ProngsPrompt > 1) {
            keepEvent[kDoubleCharm3P] = true;
          }
        }
      }
      if (n23Prongs > 1 && enableDoubleCharmChannels->get(0u, 2u)) {
        if (enableDoubleCharmChannels->get(1u, 2u)) {
          keepEvent[kDoubleCharmMix] = true;
        } else {
          if (n23ProngsPrompt > 1) {
            keepEvent[kDoubleCharmMix] = true;
          }
        }
      }

      // apply downscale factors, if required
      if (applyDownscale) {
        auto rndValue = gRandom->Rndm();
        for (int iTrigger{0}; iTrigger < kNtriggersHF; ++iTrigger) {
          if (rndValue > downscaleFactors->get(iTrigger, 0u)) {
            keepEvent[iTrigger] = false;
          }
        }
      }

      tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix], keepEvent[kV0Charm2P], keepEvent[kV0Charm3P], keepEvent[kCharmBarToXiBach], keepEvent[kSigmaCPPK], keepEvent[kSigmaC0K0], keepEvent[kPhotonCharm2P], keepEvent[kPhotonCharm3P], keepEvent[kSingleCharm2P], keepEvent[kSingleCharm3P], keepEvent[kSingleNonPromptCharm2P], keepEvent[kSingleNonPromptCharm3P], keepEvent[kCharmBarToXi2Bach], keepEvent[kPrCharm2P], keepEvent[kBtoJPsiKa], keepEvent[kBtoJPsiKstar], keepEvent[kBtoJPsiPhi], keepEvent[kBtoJPsiPrKa], keepEvent[kBtoJPsiPi]);

      if (!std::accumulate(keepEvent, keepEvent + kNtriggersHF, 0)) {
        hProcessedEvents->Fill(1);
      } else {
        for (int iTrigger{0}; iTrigger < kNtriggersHF; ++iTrigger) {
          if (keepEvent[iTrigger]) {
            hProcessedEvents->Fill(iTrigger + 2);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfFilter>(cfg));

  return workflow;
}
