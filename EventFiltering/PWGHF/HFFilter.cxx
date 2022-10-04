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
// O2 includes

/// \file HFFilter.cxx
/// \brief task for selection of events with HF signals
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "EventFiltering/filterTables.h"

#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

#include <cmath>
#include <string>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

// ML application
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace hf_cuts_single_track;
using namespace hf_cuts_bdt_multiclass;

namespace
{

enum HfTriggers {
  kHighPt2P = 0,
  kHighPt3P,
  kBeauty3P,
  kBeauty4P,
  kFemto2P,
  kFemto3P,
  kDoubleCharm2P,
  kDoubleCharm3P,
  kDoubleCharmMix,
  kNtriggersHF
};

enum charmParticles {
  kD0 = 0,
  kDplus,
  kDs,
  kLc,
  kXic,
  kNCharmParticles
};

enum beautyParticles {
  kBplus = 0,
  kB0toDStar,
  kB0,
  kBs,
  kLb,
  kXib,
  kNBeautyParticles
};

enum beautyTrackSelection {
  kRejected = 0,
  kSoftPion,
  kRegular
};

static const std::array<std::string, kNtriggersHF> HfTriggerNames{"highPt", "beauty", "femto", "doubleCharm"};
static const std::array<std::string, kNCharmParticles> charmParticleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const std::array<std::string, kNBeautyParticles> beautyParticleNames{"Bplus", "B0toDStar", "B0", "Bs", "Lb", "Xib"};

static const float massPi = RecoDecay::getMassPDG(211);
static const float massK = RecoDecay::getMassPDG(321);
static const float massProton = RecoDecay::getMassPDG(2212);
static const float massPhi = RecoDecay::getMassPDG(333);
static const float massD0 = RecoDecay::getMassPDG(421);
static const float massDPlus = RecoDecay::getMassPDG(411);
static const float massDs = RecoDecay::getMassPDG(431);
static const float massLc = RecoDecay::getMassPDG(4122);
static const float massXic = RecoDecay::getMassPDG(4232);
static const float massDStar = RecoDecay::getMassPDG(413);
static const float massBPlus = RecoDecay::getMassPDG(511);
static const float massB0 = RecoDecay::getMassPDG(521);
static const float massBs = RecoDecay::getMassPDG(531);
static const float massLb = RecoDecay::getMassPDG(5122);
static const float massXib = RecoDecay::getMassPDG(5232);

} // namespace

namespace o2::aod
{
namespace extra2Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra2Prong
namespace extra3Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra3Prong
DECLARE_SOA_TABLE(Colls2Prong, "AOD", "COLLSID2P", o2::aod::extra2Prong::CollisionId);
DECLARE_SOA_TABLE(Colls3Prong, "AOD", "COLLSID3P", o2::aod::extra3Prong::CollisionId);

namespace hftraining2p
{
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);       //!
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float); //!
DECLARE_SOA_COLUMN(PT2Prong, pT2Prong, float);         //!
DECLARE_SOA_COLUMN(PT1, pT1, float);                   //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);     //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);       //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float); //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float); //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                   //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);     //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);       //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float); //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float); //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);    //!
} // namespace hftraining2p
DECLARE_SOA_TABLE(HFTrigTrain2P, "AOD", "HFTRIGTRAIN2P", //!
                  hftraining2p::InvMassD0,
                  hftraining2p::InvMassD0bar,
                  hftraining2p::PT2Prong,
                  hftraining2p::PT1,
                  hftraining2p::DCAPrimXY1,
                  hftraining2p::DCAPrimZ1,
                  hftraining2p::NsigmaPiTPC1,
                  hftraining2p::NsigmaKaTPC1,
                  hftraining2p::NsigmaPiTOF1,
                  hftraining2p::NsigmaKaTOF1,
                  hftraining2p::PT2,
                  hftraining2p::DCAPrimXY2,
                  hftraining2p::DCAPrimZ2,
                  hftraining2p::NsigmaPiTPC2,
                  hftraining2p::NsigmaKaTPC2,
                  hftraining2p::NsigmaPiTOF2,
                  hftraining2p::NsigmaKaTOF2,
                  hftraining2p::FlagOrigin);

namespace hftraining3p
{
DECLARE_SOA_COLUMN(InvMassDplus, invMassDplus, float);           //!
DECLARE_SOA_COLUMN(InvMassDsToKKPi, invMassDsToKKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassDsToPiKK, invMassDsToPiKK, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPKPi, invMassLcToPKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPiKP, invMassLcToPiKP, float);     //!
DECLARE_SOA_COLUMN(InvMassXicToPKPi, invMassXicToPKPi, float);   //!
DECLARE_SOA_COLUMN(InvMassXicToPiKP, invMassXicToPiKP, float);   //!
DECLARE_SOA_COLUMN(PT3Prong, pT3Prong, float);                   //!
DECLARE_SOA_COLUMN(PT1, pT1, float);                             //!
DECLARE_SOA_COLUMN(DeltaMassKKFirst, deltaMassKKFirst, float);   //!
DECLARE_SOA_COLUMN(DeltaMassKKSecond, deltaMassKKSecond, float); //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC1, nsigmaPrTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF1, nsigmaPrTOF1, float);           //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC2, nsigmaPrTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF2, nsigmaPrTOF2, float);           //!
DECLARE_SOA_COLUMN(PT3, pT3, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY3, dcaPrimXY3, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ3, dcaPrimZ3, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC3, nsigmaPiTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC3, nsigmaKaTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC3, nsigmaPrTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF3, nsigmaPiTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF3, nsigmaKaTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF3, nsigmaPrTOF3, float);           //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);              //!
DECLARE_SOA_COLUMN(Channel, channel, int8_t);                    //!
DECLARE_SOA_COLUMN(HFSelBit, hfselbit, int8_t);                  //!
} // namespace hftraining3p
DECLARE_SOA_TABLE(HFTrigTrain3P, "AOD", "HFTRIGTRAIN3P", //!
                  hftraining3p::InvMassDplus,
                  hftraining3p::InvMassDsToKKPi,
                  hftraining3p::InvMassDsToPiKK,
                  hftraining3p::InvMassLcToPKPi,
                  hftraining3p::InvMassLcToPiKP,
                  hftraining3p::InvMassXicToPKPi,
                  hftraining3p::InvMassXicToPiKP,
                  hftraining3p::PT3Prong,
                  hftraining3p::DeltaMassKKFirst,
                  hftraining3p::DeltaMassKKSecond,
                  hftraining3p::PT1,
                  hftraining3p::DCAPrimXY1,
                  hftraining3p::DCAPrimZ1,
                  hftraining3p::NsigmaPiTPC1,
                  hftraining3p::NsigmaKaTPC1,
                  hftraining3p::NsigmaPrTPC1,
                  hftraining3p::NsigmaPiTOF1,
                  hftraining3p::NsigmaKaTOF1,
                  hftraining3p::NsigmaPrTOF1,
                  hftraining3p::PT2,
                  hftraining3p::DCAPrimXY2,
                  hftraining3p::DCAPrimZ2,
                  hftraining3p::NsigmaPiTPC2,
                  hftraining3p::NsigmaKaTPC2,
                  hftraining3p::NsigmaPrTPC2,
                  hftraining3p::NsigmaPiTOF2,
                  hftraining3p::NsigmaKaTOF2,
                  hftraining3p::NsigmaPrTOF2,
                  hftraining3p::PT3,
                  hftraining3p::DCAPrimXY3,
                  hftraining3p::DCAPrimZ3,
                  hftraining3p::NsigmaPiTPC3,
                  hftraining3p::NsigmaKaTPC3,
                  hftraining3p::NsigmaPrTPC3,
                  hftraining3p::NsigmaPiTOF3,
                  hftraining3p::NsigmaKaTOF3,
                  hftraining3p::NsigmaPrTOF3,
                  hftraining3p::FlagOrigin,
                  hftraining3p::Channel,
                  hftraining3p::HFSelBit);
} // namespace o2::aod

struct AddCollisionId {

  Produces<o2::aod::Colls2Prong> colls2Prong;
  Produces<o2::aod::Colls3Prong> colls3Prong;

  void process(aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::Tracks const&)
  {
    for (const auto& cand2Prong : cand2Prongs) {
      colls2Prong(cand2Prong.index0_as<aod::Tracks>().collisionId());
    }
    for (const auto& cand3Prong : cand3Prongs) {
      colls3Prong(cand3Prong.index0_as<aod::Tracks>().collisionId());
    }
  }
};

struct HfFilter { // Main struct for HF triggers

  Produces<aod::HfFilters> tags;
  Produces<aod::HFTrigTrain2P> train2P;
  Produces<aod::HFTrigTrain3P> train3P;

  Configurable<bool> activateQA{"activateQA", false, "flag to enable QA histos"};

  // parameters for high-pT triggers
  Configurable<float> pTThreshold2Prong{"pTThreshold2Prong", 8., "pT treshold for high pT 2-prong candidates for kHighPt triggers in GeV/c"};
  Configurable<float> pTThreshold3Prong{"pTThreshold3Prong", 8., "pT treshold for high pT 3-prong candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<float> deltaMassBPlus{"deltaMassBPlus", 0.3, "invariant-mass delta with respect to the B+ mass"};
  Configurable<float> deltaMassB0{"deltaMassB0", 0.3, "invariant-mass delta with respect to the B0 mass"};
  Configurable<float> deltaMassBs{"deltaMassBs", 0.3, "invariant-mass delta with respect to the Bs mass"};
  Configurable<float> deltaMassCharmHadronForBeauty{"deltaMassCharmHadronForBeauty", 0.04, "invariant-mass delta for charm"};
  Configurable<float> deltaMassLb{"deltaMassLb", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassXib{"deltaMassXib", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassDStar{"deltaMassDStar", 0.04, "invariant-mass delta with respect to the D* mass for B0 -> D*pi"};
  Configurable<float> pTMinBeautyBachelor{"pTMinBeautyBachelor", 0.5, "minumum pT for bachelor pion track used to build b-hadron candidates"};
  Configurable<float> pTMinSoftPion{"pTMinSoftPion", 0.1, "minumum pT for soft pion track used to build D* mesons in the b-hadron decay chain"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{hf_cuts_single_track::pTBinsTrack_v}, "track pT bin limits for DCAXY pT-depentend cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::npTBinsTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::pTBinLabelsTrack, hf_cuts_single_track::cutVarLabelsTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::npTBinsTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::pTBinLabelsTrack, hf_cuts_single_track::cutVarLabelsTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};
  Configurable<float> nsigmaTPCProtonLc{"nsigmaTPCProtonLc", 3., "Maximum value for TPC PID proton Nsigma for Lc"};
  Configurable<float> nsigmaTOFProtonLc{"nsigmaTOFProtonLc", 3., "Maximum value for TOF PID proton Nsigma for Lc"};

  // parameters for femto triggers
  Configurable<float> femtoMaxRelativeMomentum{"femtoMaxRelativeMomentum", 2., "Maximal allowed value for relative momentum between charm-proton pairs in GeV/c"};
  Configurable<float> femtoMinProtonPt{"femtoMinProtonPt", 0.5, "Minimal required transverse momentum for proton in GeV/c"};
  Configurable<bool> femtoProtonOnlyTOF{"femtoProtonOnlyTOF", true, "Use only TOF information for proton identification if true"};
  Configurable<float> femtoMaxNsigmaProton{"femtoMaxNsigmaProton", 3., "Maximum value for PID proton Nsigma for femto triggers"};

  // array of single-track cuts for pion
  std::array<LabeledArray<double>, 2> cutsSingleTrackBeauty;

  // parameters for production of training samples
  Configurable<bool> fillSignal{"fillSignal", true, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillBackground{"fillBackground", true, "Flag to fill derived tables with background for ML trainings"};
  Configurable<double> donwSampleBkgFactor{"donwSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};

  // parameters for ML application with ONNX
  Configurable<bool> singleThreadInference{"singleThreadInference", true, "Run ML inference single thread"};
  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::vector<double>> pTBinsBDT{"pTBinsBDT", std::vector<double>{hf_cuts_bdt_multiclass::pTBinsVec}, "track pT bin limits for BDT cut"};

  Configurable<std::string> onnxFileD0ToKPiConf{"onnxFileD0ToKPiConf", "/cvmfs/alice.cern.ch/data/analysis/2022/vAN-20220124/PWGHF/o2/trigger/XGBoostModel.onnx", "ONNX file for ML model for D0 candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreD0ToKPi{"thresholdBDTScoreD0ToKPi", {hf_cuts_bdt_multiclass::cutsBDT[0], hf_cuts_bdt_multiclass::npTBins, hf_cuts_bdt_multiclass::nCutBDTScores, hf_cuts_bdt_multiclass::pTBinLabels, hf_cuts_bdt_multiclass::cutBDTLabels}, "Threshold values for BDT output scores of D0 candidates"};

  Configurable<std::string> onnxFileDPlusToPiKPiConf{"onnxFileDPlusToPiKPiConf", "", "ONNX file for ML model for D+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreDPlusToPiKPi{"thresholdBDFScoreDPlusToPiKPi", {hf_cuts_bdt_multiclass::cutsBDT[0], hf_cuts_bdt_multiclass::npTBins, hf_cuts_bdt_multiclass::nCutBDTScores, hf_cuts_bdt_multiclass::pTBinLabels, hf_cuts_bdt_multiclass::cutBDTLabels}, "Threshold values for BDT output scores of D+ candidates"};

  Configurable<std::string> onnxFileDSToPiKKConf{"onnxFileDSToPiKKConf", "", "ONNX file for ML model for Ds+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreDSToPiKK{"thresholdBDFScoreDSToPiKK", {hf_cuts_bdt_multiclass::cutsBDT[0], hf_cuts_bdt_multiclass::npTBins, hf_cuts_bdt_multiclass::nCutBDTScores, hf_cuts_bdt_multiclass::pTBinLabels, hf_cuts_bdt_multiclass::cutBDTLabels}, "Threshold values for BDT output scores of Ds+ candidates"};

  Configurable<std::string> onnxFileLcToPiKPConf{"onnxFileLcToPiKPConf", "", "ONNX file for ML model for Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreLcToPiKP{"thresholdBDFScoreLcToPiKP", {hf_cuts_bdt_multiclass::cutsBDT[0], hf_cuts_bdt_multiclass::npTBins, hf_cuts_bdt_multiclass::nCutBDTScores, hf_cuts_bdt_multiclass::pTBinLabels, hf_cuts_bdt_multiclass::cutBDTLabels}, "Threshold values for BDT output scores of Lc+ candidates"};

  Configurable<std::string> onnxFileXicToPiKPConf{"onnxFileXicToPiKPConf", "", "ONNX file for ML model for Xic+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDFScoreXicToPiKP{"thresholdBDFScoreXicToPiKP", {hf_cuts_bdt_multiclass::cutsBDT[0], hf_cuts_bdt_multiclass::npTBins, hf_cuts_bdt_multiclass::nCutBDTScores, hf_cuts_bdt_multiclass::pTBinLabels, hf_cuts_bdt_multiclass::cutBDTLabels}, "Threshold values for BDT output scores of Xic+ candidates"};

  // array of ONNX config and BDT thresholds
  std::array<std::string, kNCharmParticles> onnxFiles;
  std::array<LabeledArray<double>, kNCharmParticles> thresholdBDTScores;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  std::shared_ptr<TH1> hProcessedEvents;

  // QA histos
  std::shared_ptr<TH1> hN2ProngCharmCand, hN3ProngCharmCand;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmHighPt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hCharmProtonKstarDistr{};
  std::array<std::shared_ptr<TH2>, kNBeautyParticles> hMassVsPtB{};
  std::array<std::shared_ptr<TH2>, kNCharmParticles + 1> hMassVsPtC{};
  std::shared_ptr<TH2> hProtonTPCPID, hProtonTOFPID;
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreBkg{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScorePrompt{};
  std::array<std::shared_ptr<TH1>, kNCharmParticles> hBDTScoreNonPrompt{};

  // ONNX
  // Ort::Env env;
  std::array<std::vector<std::string>, kNCharmParticles> inputNamesML{};
  std::array<std::vector<std::vector<int64_t>>, kNCharmParticles> inputShapesML{};
  std::array<std::vector<std::string>, kNCharmParticles> outputNamesML{};
  std::array<std::vector<std::vector<int64_t>>, kNCharmParticles> outputShapesML{};
  std::array<std::shared_ptr<Ort::Experimental::Session>, kNCharmParticles> sessionML = {nullptr, nullptr, nullptr, nullptr, nullptr};
  std::array<Ort::SessionOptions, kNCharmParticles> sessionOptions{Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions(), Ort::SessionOptions()};
  std::array<int, kNCharmParticles> dataTypeML{};
  std::array<Ort::Env, kNCharmParticles> env = {
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-d0-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-dplus-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-ds-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-lc-triggers"},
    Ort::Env{ORT_LOGGING_LEVEL_ERROR, "ml-model-xic-triggers"}};

  void init(o2::framework::InitContext&)
  {
    cutsSingleTrackBeauty = {cutsTrackBeauty3Prong, cutsTrackBeauty4Prong};

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;counts", HistType::kTH1F, {{kNtriggersHF + 2, -0.5, kNtriggersHF + 1.5}});
    std::array<std::string, kNtriggersHF + 2> eventTitles = {"all", "rejected", "w/ high-#it{p}_{T} 2p charm", "w/ high-#it{p}_{T} 3p charm", "w/ 3p beauty", "w/ 4p beauty", "w/ 2p femto", "w/ 3p femto", "w/ 2p double charm", "w/ 3p double charm", "w/ 2p and 3p double charm"};
    for (auto iBin = 0; iBin < kNtriggersHF + 2; ++iBin) {
      hProcessedEvents->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    if (activateQA) {
      hN2ProngCharmCand = registry.add<TH1>("fN2ProngCharmCand", "Number of 2-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1F, {{50, -0.5, 49.5}});
      hN3ProngCharmCand = registry.add<TH1>("fN3ProngCharmCand", "Number of 3-prong charm candidates per event;#it{N}_{candidates};counts", HistType::kTH1F, {{50, -0.5, 49.5}});
      for (int iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        hCharmHighPt[iCharmPart] = registry.add<TH1>(Form("f%sHighPt", charmParticleNames[iCharmPart].data()), Form("#it{p}_{T} distribution of triggered high-#it{p}_{T} %s candidates;#it{p}_{T} (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 50.}});
        hCharmProtonKstarDistr[iCharmPart] = registry.add<TH1>(Form("f%sProtonKstarDistr", charmParticleNames[iCharmPart].data()), Form("#it{k}* distribution of triggered p#minus%s pairs;#it{k}* (GeV/#it{c});counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
        hMassVsPtC[iCharmPart] = registry.add<TH2>(Form("fMassVsPt%s", charmParticleNames[iCharmPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", charmParticleNames[iCharmPart].data()), HistType::kTH2F, {{100, 0., 50.}, {300, 1.60, 2.60}});
        if (applyML) {
          hBDTScoreBkg[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreBkgDistr", charmParticleNames[iCharmPart].data()), Form("BDT background score distribution for %s;BDT background score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
          hBDTScorePrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScorePromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT prompt score distribution for %s;BDT prompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
          hBDTScoreNonPrompt[iCharmPart] = registry.add<TH1>(Form("f%sBDTScoreNonPromptDistr", charmParticleNames[iCharmPart].data()), Form("BDT nonprompt score distribution for %s;BDT nonprompt score;counts", charmParticleNames[iCharmPart].data()), HistType::kTH1F, {{100, 0., 1.}});
        }
      }
      hMassVsPtC[kNCharmParticles] = registry.add<TH2>("fMassVsPtDStar", "#it{M} vs. #it{p}_{T} distribution of triggered DStar candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", HistType::kTH2F, {{100, 0., 50.}, {300, 1.60, 2.60}});
      for (int iBeautyPart{0}; iBeautyPart < kNBeautyParticles; ++iBeautyPart) {
        hMassVsPtB[iBeautyPart] = registry.add<TH2>(Form("fMassVsPt%s", beautyParticleNames[iBeautyPart].data()), Form("#it{M} vs. #it{p}_{T} distribution of triggered %s candidates;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2});counts", beautyParticleNames[iBeautyPart].data()), HistType::kTH2F, {{100, 0., 50.}, {220, 4.9, 6.0}});
      }
      hProtonTPCPID = registry.add<TH2>("fProtonTPCPID", "#it{N}_{#sigma}^{TPC} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TPC}", HistType::kTH2F, {{100, 0., 10.}, {200, -10., 10.}});
      hProtonTOFPID = registry.add<TH2>("fProtonTOFPID", "#it{N}_{#sigma}^{TOF} vs. #it{p} for selected protons;#it{p} (GeV/#it{c});#it{N}_{#sigma}^{TOF}", HistType::kTH2F, {{100, 0., 10.}, {200, -10., 10.}});
    }

    // init ONNX runtime session
    if (applyML) {
      thresholdBDTScores = {
        thresholdBDTScoreD0ToKPi,
        thresholdBDFScoreDPlusToPiKPi,
        thresholdBDFScoreDSToPiKK,
        thresholdBDFScoreLcToPiKP,
        thresholdBDFScoreXicToPiKP};

      onnxFiles = {
        onnxFileD0ToKPiConf,
        onnxFileDPlusToPiKPiConf,
        onnxFileDSToPiKKConf,
        onnxFileLcToPiKPConf,
        onnxFileXicToPiKPConf};

      for (auto iCharmPart{0}; iCharmPart < kNCharmParticles; ++iCharmPart) {
        if (onnxFiles[iCharmPart] != "") {
          if (singleThreadInference) {
            sessionOptions[iCharmPart].SetIntraOpNumThreads(1);
            sessionOptions[iCharmPart].SetInterOpNumThreads(1);
          }
          sessionML[iCharmPart].reset(new Ort::Experimental::Session{env[iCharmPart], onnxFiles[iCharmPart], sessionOptions[iCharmPart]});
          inputNamesML[iCharmPart] = sessionML[iCharmPart]->GetInputNames();
          inputShapesML[iCharmPart] = sessionML[iCharmPart]->GetInputShapes();
          if (inputShapesML[iCharmPart][0][0] < 0) {
            LOGF(warning, Form("Model for %s with negative input shape likely because converted with ummingbird, setting it to 1.", charmParticleNames[iCharmPart].data()));
            inputShapesML[iCharmPart][0][0] = 1;
          }
          outputNamesML[iCharmPart] = sessionML[iCharmPart]->GetOutputNames();
          outputShapesML[iCharmPart] = sessionML[iCharmPart]->GetOutputShapes();

          Ort::TypeInfo typeInfo = sessionML[iCharmPart]->GetInputTypeInfo(0);
          auto tensorInfo = typeInfo.GetTensorTypeAndShapeInfo();
          dataTypeML[iCharmPart] = tensorInfo.GetElementType();
        }
      }
    }
  }

  /// Single-track cuts for bachelor track of beauty candidates
  /// \param track is a track
  /// \param candType candidate type (3-prong or 4-prong beauty candidate)
  /// \return 0 if track is rejected, 1 if track is soft pion, 2 if it is regular beauty
  template <typename T>
  int isSelectedTrackForBeauty(const T& track, const int candType)
  {
    auto pT = track.pt();
    auto pTBinTrack = findBin(pTBinsTrack, pT);
    if (pTBinTrack == -1) {
      return kRejected;
    }

    if (pT < pTMinSoftPion) { // soft pion should be more stringent than usual tracks
      return kRejected;
    }

    if (std::abs(track.eta()) > 0.8) {
      return kRejected;
    }

    if (std::abs(track.dcaZ()) > 2.f) {
      return kRejected;
    }

    if (std::abs(track.dcaXY()) < cutsSingleTrackBeauty[candType - 2].get(pTBinTrack, "min_dcaxytoprimary")) {
      return kRejected; // minimum DCAxy
    }
    if (std::abs(track.dcaXY()) > cutsSingleTrackBeauty[candType - 2].get(pTBinTrack, "max_dcaxytoprimary")) {
      return kRejected; // maximum DCAxy
    }

    // below only regular beauty tracks, not required for soft pions
    if (pT < pTMinBeautyBachelor) {
      return kSoftPion;
    }

    return kRegular;
  }

  /// Basic selection of proton candidates
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedProton4Femto(const T& track)
  {
    if (track.pt() < femtoMinProtonPt) {
      return false;
    }

    if (std::abs(track.eta()) > 0.8) {
      return false;
    }

    if (track.isGlobalTrack() != (uint8_t) true) {
      return false; // use only global tracks
    }

    float NSigmaTPC = track.tpcNSigmaPr();
    float NSigmaTOF = track.tofNSigmaPr();
    float NSigma;

    if (femtoProtonOnlyTOF) {
      NSigma = abs(NSigmaTOF);
    } else {
      NSigma = sqrt(NSigmaTPC * NSigmaTPC + NSigmaTOF * NSigmaTOF);
    }

    if (NSigma > femtoMaxNsigmaProton) {
      return false;
    }

    if (activateQA) {
      hProtonTPCPID->Fill(track.p(), NSigmaTPC);
      hProtonTOFPID->Fill(track.p(), NSigmaTOF);
    }

    return true;
  }

  /// Basic selection of proton candidates for Lc
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedProton4CharmBaryons(const T& track)
  {
    float NSigmaTPC = track.tpcNSigmaPr();
    float NSigmaTOF = track.tofNSigmaPr();

    if (std::abs(NSigmaTPC) > nsigmaTPCProtonLc) {
      return false;
    }
    if (track.hasTOF() && std::abs(NSigmaTOF) > nsigmaTOFProtonLc) {
      return false;
    }

    return true;
  }

  /// Basic additional selection of Ds candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \return BIT(0) for KKpi, BIT(1) for piKK
  template <typename T>
  int8_t isDsPreselected(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge)
  {
    auto invMassKKFirst = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massK, massK});
    auto invMassKKSecond = RecoDecay::m(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massK, massK});

    int8_t retValue = 0;
    if (std::abs(invMassKKFirst - massPhi) < 0.02) {
      retValue |= BIT(0);
    }
    if (std::abs(invMassKKSecond - massPhi) < 0.02) {
      retValue |= BIT(1);
    }

    return retValue;
  }

  /// Basic additional selection of Lc->pKpi and Xic->pKpi candidates
  /// \param TrackSameChargeFirst is the first same-charge track
  /// \param TrackSameChargeSecond is the second same-charge track
  /// \return BIT(0) for pKpi, BIT(1) for piKp
  template <typename T>
  int8_t isCharmBaryonPreselected(const T& TrackSameChargeFirst, const T& TrackSameChargeSecond)
  {
    int8_t retValue = 0;
    if (isSelectedProton4CharmBaryons(TrackSameChargeFirst)) {
      retValue |= BIT(0);
    }
    if (isSelectedProton4CharmBaryons(TrackSameChargeSecond)) {
      retValue |= BIT(1);
    }

    return retValue;
  }

  /// Mass selection of D0 candidates to build Bplus candidates
  /// \param pTrackPos is the positive track momentum
  /// \param pTrackNeg is the negative track momentum
  /// \param ptD is the pt of the D0 meson candidate
  /// \return 1 for D0, 2 for D0bar, 3 for both
  template <typename T>
  int8_t isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD)
  {
    auto invMassD0 = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massK});
    auto invMassD0bar = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massK, massPi});

    if (activateQA) {
      hMassVsPtC[kD0]->Fill(ptD, invMassD0);
      hMassVsPtC[kD0]->Fill(ptD, invMassD0bar);
    }

    int8_t retValue = 0;
    if (std::abs(invMassD0 - massD0) < deltaMassCharmHadronForBeauty) {
      retValue += 1;
    }
    if (std::abs(invMassD0bar - massD0) < deltaMassCharmHadronForBeauty) {
      retValue += 2;
    }

    return retValue;
  }

  /// Mass selection of D+ candidates to build B0 candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the D+ meson candidate
  /// \return BIT(0) (==1) for D+, 0 otherwise
  template <typename T>
  int8_t isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD)
  {
    auto invMassDplus = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massK});
    if (activateQA) {
      hMassVsPtC[kDplus]->Fill(ptD, invMassDplus);
    }

    if (std::abs(invMassDplus - massDPlus) > deltaMassCharmHadronForBeauty) {
      return 0;
    }

    return BIT(0);
  }

  /// Mass selection of of Ds candidates to build Bs candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the Ds meson candidate
  /// \return BIT(0) for KKpi, BIT(1) for piKK, BIT(2) for phipi, BIT(3) for piphi
  template <typename T>
  int8_t isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassDsToKKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massK, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kDs]->Fill(ptD, invMassDsToKKPi);
      }
      if (std::abs(invMassDsToKKPi - massDs) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassDsToPiKK = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massK});
      if (activateQA) {
        hMassVsPtC[kDs]->Fill(ptD, invMassDsToPiKK);
      }
      if (std::abs(invMassDsToPiKK - massDs) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Mass selection of Lc candidates to build Lb candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeSecond is the second same-charge track momentum
  /// \param pTrackOppositeCharge is the opposite charge track momentum
  /// \param ptLc is the pt of the D0 meson candidate
  /// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
  template <typename T>
  int8_t isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPKPi);
      }
      if (std::abs(invMassLcToPKPi - massLc) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
      if (activateQA) {
        hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPiKP);
      }
      if (std::abs(invMassLcToPiKP - massLc) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Mass selection of Xic candidates to build Lb candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeSecond is the second same-charge track momentum
  /// \param pTrackOppositeCharge is the opposite charge track momentum
  /// \param ptXic is the pt of the D0 meson candidate
  /// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
  template <typename T>
  int8_t isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const int8_t isSelected)
  {
    int8_t retValue = 0;
    if (TESTBIT(isSelected, 0)) {
      auto invMassXicToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
      if (activateQA) {
        hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPKPi);
      }
      if (std::abs(invMassXicToPKPi - massXic) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(0);
      }
    }
    if (TESTBIT(isSelected, 1)) {
      auto invMassXicToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});
      if (activateQA) {
        hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPiKP);
      }
      if (std::abs(invMassXicToPiKP - massXic) < deltaMassCharmHadronForBeauty) {
        retValue |= BIT(1);
      }
    }

    return retValue;
  }

  /// Single-track cuts for bachelor track of beauty candidates
  /// \param scores is a 3-element array with BDT out scores
  /// \return 0 if rejected, otherwise bitmap with BIT(RecoDecay::OriginType::Prompt) and/or BIT(RecoDecay::OriginType::NonPrompt) on
  template <typename T>
  int8_t isBDTSelected(const T& scores, const int candType)
  {
    int8_t retValue = 0;

    if (scores[0] > thresholdBDTScores[candType].get(0u, "BDTbkg")) {
      return retValue;
    }
    if (scores[1] < thresholdBDTScores[candType].get(0u, "BDTprompt")) {
      retValue |= BIT(RecoDecay::OriginType::Prompt);
    }
    if (scores[2] < thresholdBDTScores[candType].get(0u, "BDTnonprompt")) {
      retValue |= BIT(RecoDecay::OriginType::NonPrompt);
    }

    return retValue;
  }

  /// Computation of the relative momentum between particle pairs
  /// \param track is a track
  /// \param ProtonMass is the mass of a proton
  /// \param CharmCandMomentum is the three momentum of a charm candidate
  /// \param CharmMass is the mass of the charm hadron
  /// \return relative momentum of pair
  template <typename T> // template <typename T, typename C>
  float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)
  {
    ROOT::Math::PxPyPzMVector part1(track.px(), track.py(), track.pz(), massProton);
    ROOT::Math::PxPyPzMVector part2(CharmCandMomentum[0], CharmCandMomentum[1], CharmCandMomentum[2], CharmMass);

    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);
    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;

    float kStar = 0.5 * trackRelK.P();
    return kStar;
  } // float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

  using HfTrackIndexProng2withColl = soa::Join<aod::Hf2Prongs, aod::Colls2Prong>;
  using HfTrackIndexProng3withColl = soa::Join<aod::Hf3Prongs, aod::Colls3Prong>;
  using BigTracksMCPID = soa::Join<aod::BigTracksExtended, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::BigTracksMC>;

  Filter trackFilter = requireGlobalTrackWoDCAInFilter();
  using BigTracksWithProtonPID = soa::Filtered<soa::Join<aod::BigTracksExtended, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  void process(aod::Collision const& collision,
               HfTrackIndexProng2withColl const& cand2Prongs,
               HfTrackIndexProng3withColl const& cand3Prongs,
               BigTracksWithProtonPID const& tracks)
  {
    hProcessedEvents->Fill(0);

    // collision process loop
    bool keepEvent[kNtriggersHF]{false};
    //

    int n2Prongs{0}, n3Prongs{0};

    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      if (!TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_prong2::DecayType::D0ToPiK)) { // check if it's a D0
        continue;
      }

      auto trackPos = cand2Prong.index0_as<BigTracksWithProtonPID>(); // positive daughter
      auto trackNeg = cand2Prong.index1_as<BigTracksWithProtonPID>(); // negative daughter
      std::array<float, 3> pVecPos = {trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};

      bool isCharmTagged{true}, isBeautyTagged{true};

      // apply ML models
      if (applyML && onnxFiles[kD0] != "") {
        isCharmTagged = false;
        isBeautyTagged = false;
        // TODO: add more feature configurations
        std::vector<Ort::Value> inputTensorD0;
        std::vector<float> inputFeaturesD0{trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::vector<double> inputFeaturesDoD0{trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ()};
        if (dataTypeML[kD0] == 1) {
          inputTensorD0.push_back(Ort::Experimental::Value::CreateTensor<float>(inputFeaturesD0.data(), inputFeaturesD0.size(), inputShapesML[kD0][0]));
        } else if (dataTypeML[kD0] == 11) {
          inputTensorD0.push_back(Ort::Experimental::Value::CreateTensor<double>(inputFeaturesDoD0.data(), inputFeaturesDoD0.size(), inputShapesML[kD0][0]));
        } else {
          LOG(fatal) << "Error running model inference: Unexpected input data type.";
        }

        // double-check the dimensions of the input tensor
        if (inputTensorD0[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
          assert(inputTensorD0[0].IsTensor() && inputTensorD0[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapesML[kD0][0]);
        }
        try {
          auto outputTensorD0 = sessionML[kD0]->Run(inputNamesML[kD0], inputTensorD0, outputNamesML[kD0]);
          assert(outputTensorD0.size() == outputNamesML[kD0].size() && outputTensorD0[1].IsTensor());
          auto typeInfo = outputTensorD0[1].GetTensorTypeAndShapeInfo();
          assert(typeInfo.GetElementCount() == 3); // we need multiclass
          auto scores = outputTensorD0[1].GetTensorMutableData<float>();

          if (applyML && activateQA) {
            hBDTScoreBkg[kD0]->Fill(scores[0]);
            hBDTScorePrompt[kD0]->Fill(scores[1]);
            hBDTScoreNonPrompt[kD0]->Fill(scores[2]);
          }

          int tagBDT = isBDTSelected(scores, kD0);
          isCharmTagged = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
          isBeautyTagged = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);
        } catch (const Ort::Exception& exception) {
          //LOG(error) << "Error running model inference: " << exception.what();
        }
      }

      if (!isCharmTagged && !isBeautyTagged) {
        continue;
      }

      auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::pt(pVec2Prong);

      auto selD0 = isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong);

      if (pt2Prong >= pTThreshold2Prong) {
        keepEvent[kHighPt2P] = true;
        if (activateQA) {
          hCharmHighPt[kD0]->Fill(pt2Prong);
        }
      } // end high-pT selection

      if (isCharmTagged) {
        n2Prongs++;
      } // end multi-charm selection

      for (const auto& track : tracks) { // start loop over tracks
        if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecThird = {track.px(), track.py(), track.pz()};

        if (!keepEvent[kBeauty3P] && isBeautyTagged) {
          int isTrackSelected = isSelectedTrackForBeauty(track, kBeauty3P);
          if (isTrackSelected && (((selD0 == 1 || selD0 == 3) && track.signed1Pt() < 0) || (selD0 >= 2 && track.signed1Pt() > 0))) {
            auto massCand = RecoDecay::m(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
            auto pVecBeauty3Prong = RecoDecay::pVec(pVec2Prong, pVecThird);
            auto ptCand = RecoDecay::pt(pVecBeauty3Prong);
            if (isTrackSelected == kRegular && std::abs(massCand - massBPlus) <= deltaMassBPlus) {
              keepEvent[kBeauty3P] = true;
              if (activateQA) {
                hMassVsPtB[kBplus]->Fill(ptCand, massCand);
              }
            } else if (std::abs(massCand - massDStar) <= deltaMassDStar) { // additional check for B0->D*pi polarization studies
              if (activateQA) {
                hMassVsPtC[kNCharmParticles]->Fill(ptCand, massCand);
              }
              for (const auto& trackB : tracks) { // start loop over tracks
                if (track.signed1Pt() * trackB.signed1Pt() < 0 && isSelectedTrackForBeauty(trackB, kBeauty3P) == kRegular) {
                  std::array<float, 3> pVecFourth = {trackB.px(), trackB.py(), trackB.pz()};
                  auto massCandB0 = RecoDecay::m(std::array{pVec2Prong, pVecThird, pVecFourth}, std::array{massD0, massPi, massPi});
                  if (std::abs(massCandB0 - massB0) <= deltaMassB0) {
                    keepEvent[kBeauty3P] = true;
                    if (activateQA) {
                      auto pVecBeauty4Prong = RecoDecay::pVec(pVec2Prong, pVecThird, pVecFourth);
                      auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                      hMassVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, massCandB0);
                    }
                  }
                }
              }
            }
          }
        } // end beauty selection

        // 2-prong femto
        if (!keepEvent[kFemto2P] && isCharmTagged) {
          bool isProton = isSelectedProton4Femto(track);
          if (isProton) {
            float relativeMomentum = computeRelativeMomentum(track, pVec2Prong, massD0);
            if (relativeMomentum < femtoMaxRelativeMomentum) {
              keepEvent[kFemto2P] = true;
              if (activateQA) {
                hCharmProtonKstarDistr[kD0]->Fill(relativeMomentum);
              }
            }
          }
        } // end femto selection

      } // end loop over tracks
    }   // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      std::array<int8_t, kNCharmParticles - 1> is3Prong = {
        TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DPlusToPiKPi),
        TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DsToKKPi),
        TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::LcToPKPi),
        TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::XicToPKPi)};
      if (!std::accumulate(is3Prong.begin(), is3Prong.end(), 0)) { // check if it's a D+, Ds+, Lc+ or Xic+
        continue;
      }

      auto trackFirst = cand3Prong.index0_as<BigTracksWithProtonPID>();
      auto trackSecond = cand3Prong.index1_as<BigTracksWithProtonPID>();
      auto trackThird = cand3Prong.index2_as<BigTracksWithProtonPID>();

      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      if (is3Prong[1]) { // Ds preselections
        is3Prong[1] = isDsPreselected(pVecFirst, pVecThird, pVecSecond);
      }
      if (is3Prong[2] || is3Prong[3]) { // charm baryon preselections
        auto presel = isCharmBaryonPreselected(trackFirst, trackThird);
        if (is3Prong[2]) {
          is3Prong[2] = presel;
        }
        if (is3Prong[3]) {
          is3Prong[3] = presel;
        }
      }

      std::array<int8_t, kNCharmParticles - 1> isCharmTagged = is3Prong;
      std::array<int8_t, kNCharmParticles - 1> isBeautyTagged = is3Prong;

      // apply ML models
      if (applyML) {
        isCharmTagged = std::array<int8_t, kNCharmParticles - 1>{0};
        isBeautyTagged = std::array<int8_t, kNCharmParticles - 1>{0};

        // TODO: add more feature configurations
        std::vector<float> inputFeatures{trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ()};
        std::vector<double> inputFeaturesD{trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ()};
        for (auto iCharmPart{0}; iCharmPart < kNCharmParticles - 1; ++iCharmPart) {
          if (!is3Prong[iCharmPart] || onnxFiles[iCharmPart + 1] == "") {
            continue;
          }

          std::vector<Ort::Value> inputTensor;
          if (dataTypeML[iCharmPart + 1] == 1) {
            inputTensor.push_back(Ort::Experimental::Value::CreateTensor<float>(inputFeatures.data(), inputFeatures.size(), inputShapesML[iCharmPart + 1][0]));
          } else if (dataTypeML[iCharmPart + 1] == 11) {
            inputTensor.push_back(Ort::Experimental::Value::CreateTensor<double>(inputFeaturesD.data(), inputFeaturesD.size(), inputShapesML[iCharmPart + 1][0]));
          } else {
            LOG(error) << "Error running model inference: Unexpected input data type.";
          }

          // double-check the dimensions of the input tensor
          if (inputTensor[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
            assert(inputTensor[0].IsTensor() && inputTensor[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapesML[iCharmPart + 1][0]);
          }
          try {
            auto outputTensor = sessionML[iCharmPart + 1]->Run(inputNamesML[iCharmPart + 1], inputTensor, outputNamesML[iCharmPart + 1]);
            assert(outputTensor.size() == outputNamesML[iCharmPart + 1].size() && outputTensor[1].IsTensor());
            auto typeInfo = outputTensor[1].GetTensorTypeAndShapeInfo();
            assert(typeInfo.GetElementCount() == 3); // we need multiclass
            auto scores = outputTensor[1].GetTensorMutableData<float>();

            if (applyML && activateQA) {
              hBDTScoreBkg[iCharmPart + 1]->Fill(scores[0]);
              hBDTScorePrompt[iCharmPart + 1]->Fill(scores[1]);
              hBDTScoreNonPrompt[iCharmPart + 1]->Fill(scores[2]);
            }

            int tagBDT = isBDTSelected(scores, iCharmPart + 1);
            isCharmTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::Prompt);
            isBeautyTagged[iCharmPart] = TESTBIT(tagBDT, RecoDecay::OriginType::NonPrompt);
          } catch (const Ort::Exception& exception) {
            // LOG(error) << "Error running model inference: " << exception.what();
          }
        }
      }

      if (!std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0) && !std::accumulate(isBeautyTagged.begin(), isBeautyTagged.end(), 0)) {
        continue;
      }

      if (std::accumulate(isCharmTagged.begin(), isCharmTagged.end(), 0)) {
        n3Prongs++;
      } // end multiple 3-prong selection

      auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::pt(pVec3Prong);
      float sign3Prong = trackFirst.signed1Pt() * trackSecond.signed1Pt() * trackThird.signed1Pt();

      std::array<int8_t, kNCharmParticles - 1> is3ProngInMass{0};
      if (is3Prong[0]) {
        is3ProngInMass[0] = isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
      }
      if (is3Prong[1]) {
        is3ProngInMass[1] = isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[1]);
      }
      if (is3Prong[2]) {
        is3ProngInMass[2] = isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[2]);
      }
      if (is3Prong[3]) {
        is3ProngInMass[3] = isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong, is3Prong[3]);
      }

      if (pt3Prong >= pTThreshold3Prong) {
        keepEvent[kHighPt3P] = true;
        if (activateQA) {
          for (auto iCharmPart{1}; iCharmPart < kNCharmParticles; ++iCharmPart) {
            if (is3Prong[iCharmPart - 1] && (isCharmTagged[iCharmPart - 1] || isBeautyTagged[iCharmPart - 1])) {
              hCharmHighPt[iCharmPart]->Fill(pt3Prong);
            }
          }
        }
      } // end high-pT selection

      for (const auto& track : tracks) { // start loop over tracks

        if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecFourth = {track.px(), track.py(), track.pz()};

        float massCharmHypos[kNBeautyParticles - 2] = {massDPlus, massDs, massLc, massXic};
        float massBeautyHypos[kNBeautyParticles - 2] = {massB0, massBs, massLb, massXib};
        float deltaMassHypos[kNBeautyParticles - 2] = {deltaMassB0, deltaMassBs, deltaMassLb, deltaMassXib};
        if (track.signed1Pt() * sign3Prong < 0 && isSelectedTrackForBeauty(track, kBeauty4P) == kRegular) {
          for (int iHypo{0}; iHypo < kNBeautyParticles - 2 && !keepEvent[kBeauty4P]; ++iHypo) {
            if (isBeautyTagged[iHypo] && (TESTBIT(is3ProngInMass[iHypo], 0) || TESTBIT(is3ProngInMass[iHypo], 1))) {
              auto massCandB = RecoDecay::m(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
              if (std::abs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
                keepEvent[kBeauty4P] = true;
                if (activateQA) {
                  auto pVecBeauty4Prong = RecoDecay::pVec(pVec3Prong, pVecFourth);
                  auto ptCandBeauty4Prong = RecoDecay::pt(pVecBeauty4Prong);
                  hMassVsPtB[iHypo + 2]->Fill(ptCandBeauty4Prong, massCandB);
                }
              }
            }
          }
        } // end beauty selection

        // 3-prong femto
        if (isSelectedProton4Femto(track)) {
          for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto3P]; ++iHypo) {
            if (isCharmTagged[iHypo]) {
              float relativeMomentum = computeRelativeMomentum(track, pVec3Prong, massCharmHypos[iHypo]);
              if (relativeMomentum < femtoMaxRelativeMomentum) {
                keepEvent[kFemto3P] = true;
                if (activateQA) {
                  hCharmProtonKstarDistr[iHypo + 1]->Fill(relativeMomentum);
                }
              }
            }
          }
        } // end femto selection

      } // end loop over tracks
    }   // end loop over 3-prong candidates

    if (activateQA) {
      hN2ProngCharmCand->Fill(n2Prongs);
      hN3ProngCharmCand->Fill(n3Prongs);
    }

    if (n2Prongs > 1) {
      keepEvent[kDoubleCharm2P] = true;
    }
    if (n3Prongs > 1) {
      keepEvent[kDoubleCharm3P] = true;
    }
    if (n2Prongs > 0 && n3Prongs > 0) {
      keepEvent[kDoubleCharmMix] = true;
    }

    tags(keepEvent[kHighPt2P], keepEvent[kHighPt3P], keepEvent[kBeauty3P], keepEvent[kBeauty4P], keepEvent[kFemto2P], keepEvent[kFemto3P], keepEvent[kDoubleCharm2P], keepEvent[kDoubleCharm3P], keepEvent[kDoubleCharmMix]);

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

  void
    processTraining(aod::Hf2Prongs const& cand2Prongs,
                    aod::Hf3Prongs const& cand3Prongs,
                    aod::McParticles const& particlesMC,
                    BigTracksMCPID const&)
  {
    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      auto trackPos = cand2Prong.index0_as<BigTracksMCPID>(); // positive daughter
      auto trackNeg = cand2Prong.index1_as<BigTracksMCPID>(); // negative daughter

      std::array<float, 3> pVecPos = {trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};
      auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::pt(pVec2Prong);

      auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massK});
      auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massK, massPi});

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;

      // D0(bar) → π± K∓
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, std::array{trackPos, trackNeg}, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      double pseudoRndm = trackPos.pt() * 1000. - (long)(trackPos.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train2P(invMassD0, invMassD0bar, pt2Prong, trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), trackPos.tofNSigmaPi(), trackPos.tofNSigmaKa(),
                trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ(), trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaPi(), trackNeg.tofNSigmaKa(),
                flag);
      }
    } // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      auto trackFirst = cand3Prong.index0_as<BigTracksMCPID>();  // first daughter
      auto trackSecond = cand3Prong.index1_as<BigTracksMCPID>(); // second daughter
      auto trackThird = cand3Prong.index2_as<BigTracksMCPID>();  // third daughter
      auto arrayDaughters = std::array{trackFirst, trackSecond, trackThird};

      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::pt(pVec3Prong);

      auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massPi});

      auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massK, massK, massPi});
      auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massK});

      auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      auto invMassXicToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassXicToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      float deltaMassKKFirst = -1.f;
      float deltaMassKKSecond = -1.f;
      if (TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DsToKKPi)) {
        deltaMassKKFirst = std::abs(RecoDecay::m(std::array{pVecFirst, pVecSecond}, std::array{massK, massK}) - massPhi);
        deltaMassKKSecond = std::abs(RecoDecay::m(std::array{pVecThird, pVecSecond}, std::array{massK, massK}) - massPhi);
      }
      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      int8_t channel = -1;

      // D± → π± K∓ π±
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec >= 0) {
        channel = kDplus;
      }
      if (indexRec < 0) {
        // Ds± → K± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 431, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2); // TODO: replace hard coded pdg code
        if (indexRec >= 0) {
          channel = kDs;
        }
      }
      if (indexRec < 0) {
        // Λc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kLc;
        }
      }
      if (indexRec < 0) {
        // Ξc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kXic;
        }
      }

      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      double pseudoRndm = trackFirst.pt() * 1000. - (long)(trackFirst.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train3P(invMassDplus, invMassDsToKKPi, invMassDsToPiKK, invMassLcToPKPi, invMassLcToPiKP, invMassXicToPKPi, invMassXicToPiKP, pt3Prong, deltaMassKKFirst, deltaMassKKSecond, trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackFirst.tpcNSigmaPi(), trackFirst.tpcNSigmaKa(), trackFirst.tpcNSigmaPr(), trackFirst.tofNSigmaPi(), trackFirst.tofNSigmaKa(), trackFirst.tofNSigmaPr(),
                trackSecond.pt(), trackSecond.dcaXY(), trackSecond.dcaZ(), trackSecond.tpcNSigmaPi(), trackSecond.tpcNSigmaKa(), trackSecond.tpcNSigmaPr(), trackSecond.tofNSigmaPi(), trackSecond.tofNSigmaKa(), trackSecond.tofNSigmaPr(),
                trackThird.pt(), trackThird.dcaXY(), trackThird.dcaZ(), trackThird.tpcNSigmaPi(), trackThird.tpcNSigmaKa(), trackThird.tpcNSigmaPr(), trackThird.tofNSigmaPi(), trackThird.tofNSigmaKa(), trackThird.tofNSigmaPr(),
                flag, channel, cand3Prong.hfflag());
      }
    } // end loop over 3-prong candidates
  }

  PROCESS_SWITCH(HfFilter, processTraining, "Process MC for training samples", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<AddCollisionId>(cfg));
  workflow.push_back(adaptAnalysisTask<HfFilter>(cfg));

  return workflow;
}
