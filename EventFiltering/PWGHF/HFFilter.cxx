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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "../filterTables.h"

#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

#include <cmath>
#include <string>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

// ML application
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace hf_cuts_single_track;

namespace
{

enum HfTriggers {
  kHighPt = 0,
  kBeauty,
  kFemto,
  kDoubleCharm,
  kNtriggersHF
};

enum BeautyCandType {
  kBeauty3Prong = 0, // combination of charm 2-prong and pion
  kBeauty4Prong      // combination of charm 3-prong and pion
};

enum candOrig {
  kPrompt = 0,
  kNonPrompt,
  kBkg
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
DECLARE_SOA_COLUMN(PT1, pT1, float);                   //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);     //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);       //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTPC1, nsigmaPrTPC1, float); //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTOF1, nsigmaPrTOF1, float); //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                   //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);     //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);       //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTPC2, nsigmaPrTPC2, float); //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTOF2, nsigmaPrTOF2, float); //!
DECLARE_SOA_COLUMN(PT3, pT3, float);                   //!
DECLARE_SOA_COLUMN(DCAPrimXY3, dcaPrimXY3, float);     //!
DECLARE_SOA_COLUMN(DCAPrimZ3, dcaPrimZ3, float);       //!
DECLARE_SOA_COLUMN(NsigmaPiTPC3, nsigmaPiTPC3, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTPC3, nsigmaKaTPC3, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTPC3, nsigmaPrTPC3, float); //!
DECLARE_SOA_COLUMN(NsigmaPiTOF3, nsigmaPiTOF3, float); //!
DECLARE_SOA_COLUMN(NsigmaKaTOF3, nsigmaKaTOF3, float); //!
DECLARE_SOA_COLUMN(NsigmaPrTOF3, nsigmaPrTOF3, float); //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);    //!
DECLARE_SOA_COLUMN(Channel, channel, int8_t);          //!
DECLARE_SOA_COLUMN(HFSelBit, hfselbit, int8_t);        //!
} // namespace hftraining3p
DECLARE_SOA_TABLE(HFTrigTrain3P, "AOD", "HFTRIGTRAIN3P", //!
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

  void process(aod::Hf2Prong const& cand2Prongs,
               aod::Hf3Prong const& cand3Prongs,
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
  Configurable<float> pTThreshold2Prong{"pTThreshold2Prong", 5., "pT treshold for high pT 2-prong candidates for kHighPt triggers in GeV/c"};
  Configurable<float> pTThreshold3Prong{"pTThreshold3Prong", 5., "pT treshold for high pT 3-prong candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<float> deltaMassBPlus{"deltaMassBPlus", 0.3, "invariant-mass delta with respect to the B+ mass"};
  Configurable<float> deltaMassB0{"deltaMassB0", 0.3, "invariant-mass delta with respect to the B0 mass"};
  Configurable<float> deltaMassBs{"deltaMassBs", 0.3, "invariant-mass delta with respect to the Bs mass"};
  Configurable<float> deltaMassLb{"deltaMassLb", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassXib{"deltaMassXib", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> deltaMassDStar{"deltaMassDStar", 0.1, "invariant-mass delta with respect to the D* mass for B0 -> D*pi"};
  Configurable<float> pTMinBeautyBachelor{"pTMinBeautyBachelor", 0.5, "minumum pT for bachelor pion track used to build b-hadron candidates"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{pTBinsTrack_v}, "track pT bin limits for DCAXY pT-depentend cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {cutsTrack[0], npTBinsTrack, nCutVarsTrack, pTBinLabelsTrack, cutVarLabelsTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {cutsTrack[0], npTBinsTrack, nCutVarsTrack, pTBinLabelsTrack, cutVarLabelsTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};

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
  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::string> onnxFile2ProngConf{"onnxFile2ProngConf", "/cvmfs/alice.cern.ch/data/analysis/2022/vAN-20220124/PWGHF/o2/trigger/XGBoostModel_D0ToKPi.onnx", "ONNX file for ML model for charm 2-prong candidates"};
  Configurable<float> thresholdBkgScore2Prong{"thresholdBkgScore2Prong", 0.4, "Threshold value for BDT output score on background candidates"};
  Configurable<float> thresholdPromptScore2Prong{"thresholdPromptScore2Prong", 0.5, "Threshold value for BDT output score on prompt candidates"};
  Configurable<float> thresholdNonpromptScore2Prong{"thresholdNonpromptScore2Prong", 0.1, "Threshold value for BDT output score on nonprompt candidates"};

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
  std::vector<std::string> inputNamesML2P{};
  std::vector<std::vector<int64_t>> inputShapesML2P{};
  std::vector<std::string> outputNamesML2P{};
  std::vector<std::vector<int64_t>> outputShapesML2P{};
  std::shared_ptr<Ort::Experimental::Session> session2Prong = nullptr;
  Ort::SessionOptions sessionOptions;
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "ml-model-hf-triggers"};

  void init(o2::framework::InitContext&)
  {
    cutsSingleTrackBeauty = {cutsTrackBeauty3Prong, cutsTrackBeauty4Prong};

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;counts", HistType::kTH1F, {{6, -0.5, 5.5}});
    std::array<std::string, 6> eventTitles = {"all", "rejected", "w/ high-#it{p}_{T} candidate", "w/ beauty candidate", "w/ femto candidate", "w/ double charm"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
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
      std::string onnxFile2Prong = (std::string)onnxFile2ProngConf;
      session2Prong.reset(new Ort::Experimental::Session{env, onnxFile2Prong, sessionOptions});

      inputNamesML2P = session2Prong->GetInputNames();
      inputShapesML2P = session2Prong->GetInputShapes();
      if (inputShapesML2P[0][0] < 0) {
        LOG(warning) << "Model with negative input shape likely because converted with ummingbird, setting it to 1.";
        inputShapesML2P[0][0] = 1;
      }
      outputNamesML2P = session2Prong->GetOutputNames();
      outputShapesML2P = session2Prong->GetOutputShapes();
    }
  }

  /// Single-track cuts for bachelor track of beauty candidates
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedTrackForBeauty(const T& track, const int candType)
  {
    auto pT = track.pt();
    if (pT < pTMinBeautyBachelor) {
      return false;
    }
    auto pTBinTrack = findBin(pTBinsTrack, pT);
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(track.eta()) > 0.8) {
      return false;
    }

    if (track.isGlobalTrack() != (uint8_t) true) {
      return false; // use only global tracks
    }

    if (std::abs(track.dcaXY()) < cutsSingleTrackBeauty[candType].get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; // minimum DCAxy
    }
    if (std::abs(track.dcaXY()) > cutsSingleTrackBeauty[candType].get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; // maximum DCAxy
    }
    return true;
  }

  /// Basic selection of proton candidates
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedProton(const T& track)
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

  /// Basic additional selection of D0 candidates
  /// \param pTrackPos is the positive track momentum
  /// \param pTrackNeg is the negative track momentum
  /// \param ptD is the pt of the D0 meson candidate
  /// \return 1 for D0, 2 for D0bar, 3 for both
  template <typename T>
  int isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD)
  {
    auto invMassD0 = RecoDecay::M(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massK});
    auto invMassD0bar = RecoDecay::M(std::array{pTrackPos, pTrackNeg}, std::array{massK, massPi});

    if (activateQA) {
      hMassVsPtC[kD0]->Fill(ptD, invMassD0);
      hMassVsPtC[kD0]->Fill(ptD, invMassD0bar);
    }

    int retValue = 0;
    if (std::abs(invMassD0 - massD0) < 0.04) {
      retValue += 1;
    }
    if (std::abs(invMassD0bar - massD0) < 0.04) {
      retValue += 2;
    }

    return retValue;
  }

  /// Basic additional selection of D+ candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the D+ meson candidate
  /// \return true for D+
  template <typename T>
  bool isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD)
  {
    auto invMassDplus = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massK});
    if (activateQA) {
      hMassVsPtC[kDplus]->Fill(ptD, invMassDplus);
    }

    if (std::abs(invMassDplus - massDPlus) > 0.04) {
      return false;
    }

    return true;
  }

  /// Basic additional selection of Ds candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the Ds meson candidate
  /// \return BIT(0) for KKpi, BIT(1) for piKK, BIT(2) for phipi, BIT(3) for piphi
  template <typename T>
  int isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD)
  {
    auto invMassKKFirst = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massK, massK});
    auto invMassKKSecond = RecoDecay::M(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massK, massK});

    auto invMassDsToKKPi = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massK, massK, massPi});
    auto invMassDsToPiKK = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massK});

    if (activateQA) {
      hMassVsPtC[kDs]->Fill(ptD, invMassDsToKKPi);
      hMassVsPtC[kDs]->Fill(ptD, invMassDsToPiKK);
    }

    int retValue = 0;
    if (std::abs(invMassDsToKKPi - massDs) < 0.04) {
      retValue |= BIT(0);
    }
    if (std::abs(invMassDsToPiKK - massDs) < 0.04) {
      retValue |= BIT(1);
    }
    if (std::abs(invMassKKFirst - massPhi) < 0.02) {
      retValue |= BIT(2);
    }
    if (std::abs(invMassKKSecond - massPhi) < 0.02) {
      retValue |= BIT(3);
    }

    return retValue;
  }

  /// Basic additional selection of Lc candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptLc is the pt of the D0 meson candidate
  /// \return 1 for pKpi, 2 for piKp, 3 for both
  template <typename T>
  int isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc)
  {
    auto invMassLcToPKPi = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
    auto invMassLcToPiKP = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});

    if (activateQA) {
      hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPKPi);
      hMassVsPtC[kLc]->Fill(ptLc, invMassLcToPiKP);
    }

    int retValue = 0;
    if (std::abs(invMassLcToPKPi - massLc) < 0.04) {
      retValue += 1;
    }
    if (std::abs(invMassLcToPiKP - massLc) < 0.04) {
      retValue += 2;
    }

    return retValue;
  }

  /// Basic additional selection of Xic candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptXic is the pt of the D0 meson candidate
  /// \return 1 for pKpi, 2 for piKp, 3 for both
  template <typename T>
  int isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic)
  {
    auto invMassXicToPKPi = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massK, massPi});
    auto invMassXicToPiKP = RecoDecay::M(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massK, massProton});

    if (activateQA) {
      hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPKPi);
      hMassVsPtC[kXic]->Fill(ptXic, invMassXicToPiKP);
    }

    int retValue = 0;
    if (std::abs(invMassXicToPKPi - massXic) < 0.04) {
      retValue += 1;
    }
    if (std::abs(invMassXicToPiKP - massXic) < 0.04) {
      retValue += 2;
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

  using HfTrackIndexProng2withColl = soa::Join<aod::Hf2Prong, aod::Colls2Prong>;
  using HfTrackIndexProng3withColl = soa::Join<aod::Hf3Prong, aod::Colls3Prong>;
  using BigTracksWithProtonPID = soa::Join<aod::BigTracksExtended, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using BigTracksMCPID = soa::Join<aod::BigTracksExtended, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::BigTracksMC>;

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
      if (applyML) {
        // TODO: add more feature configurations
        std::vector<float> inputFeatures2P{trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ()};
        std::vector<Ort::Value> inputTensor2P;
        inputTensor2P.push_back(Ort::Experimental::Value::CreateTensor<float>(inputFeatures2P.data(), inputFeatures2P.size(), inputShapesML2P[0]));

        // double-check the dimensions of the input tensor
        if (inputTensor2P[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
          assert(inputTensor2P[0].IsTensor() && inputTensor2P[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapesML2P[0]);
        }
        try {
          auto outputTensor2P = session2Prong->Run(inputNamesML2P, inputTensor2P, outputNamesML2P);
          assert(outputTensor2P.size() == outputNamesML2P.size() && outputTensor2P[1].IsTensor());
          auto typeInfo = outputTensor2P[1].GetTensorTypeAndShapeInfo();
          assert(typeInfo.GetElementCount() == 3); // we need multiclass
          auto scores = outputTensor2P[1].GetTensorMutableData<float>();

          if (applyML) {
            hBDTScoreBkg[kD0]->Fill(scores[0]);
            hBDTScorePrompt[kD0]->Fill(scores[1]);
            hBDTScoreNonPrompt[kD0]->Fill(scores[2]);
          }

          if (scores[0] > thresholdBkgScore2Prong) {
            continue;
          }
          if (scores[1] < thresholdPromptScore2Prong) {
            isCharmTagged = false;
          }
          if (scores[2] < thresholdNonpromptScore2Prong) {
            isBeautyTagged = false;
          }
        } catch (const Ort::Exception& exception) {
          LOG(error) << "Error running model inference: " << exception.what();
        }
      }

      if (!isCharmTagged && !isBeautyTagged) {
        continue;
      }

      auto pVec2Prong = RecoDecay::PVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::Pt(pVec2Prong);

      auto selD0 = isSelectedD0InMassRange(pVecPos, pVecNeg, pt2Prong);

      if (pt2Prong >= pTThreshold2Prong) {
        keepEvent[kHighPt] = true;
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

        if (!keepEvent[kBeauty] && isBeautyTagged) {
          if (isSelectedTrackForBeauty(track, kBeauty3Prong) && (((selD0 == 1 || selD0 == 3) && track.signed1Pt() < 0) || (selD0 >= 2 && track.signed1Pt() > 0))) {
            auto massCand = RecoDecay::M(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi});
            auto pVecBeauty3Prong = RecoDecay::PVec(pVec2Prong, pVecThird);
            auto ptCand = RecoDecay::Pt(pVecBeauty3Prong);
            if (std::abs(massCand - massBPlus) <= deltaMassBPlus) {
              keepEvent[kBeauty] = true;
              if (activateQA) {
                hMassVsPtB[kBplus]->Fill(ptCand, massCand);
              }
            } else if (std::abs(massCand - massDStar) <= deltaMassDStar) { // additional check for B0->D*pi polarization studies
              if (activateQA) {
                hMassVsPtC[kNCharmParticles]->Fill(ptCand, massCand);
              }
              for (const auto& trackB : tracks) { // start loop over tracks
                if (track.signed1Pt() * trackB.signed1Pt() < 0 && isSelectedTrackForBeauty(trackB, kBeauty3Prong)) {
                  std::array<float, 3> pVecFourth = {trackB.px(), trackB.py(), trackB.pz()};
                  auto massCandB0 = RecoDecay::M(std::array{pVec2Prong, pVecThird, pVecFourth}, std::array{massD0, massPi, massPi});
                  if (std::abs(massCandB0 - massB0) <= deltaMassB0) {
                    keepEvent[kBeauty] = true;
                    if (activateQA) {
                      auto pVecBeauty4Prong = RecoDecay::PVec(pVec2Prong, pVecThird, pVecFourth);
                      auto ptCandBeauty4Prong = RecoDecay::Pt(pVecBeauty4Prong);
                      hMassVsPtB[kB0toDStar]->Fill(ptCandBeauty4Prong, massCandB0);
                    }
                  }
                }
              }
            }
          }
        } // end beauty selection

        // 2-prong femto
        if (!keepEvent[kFemto] && isCharmTagged) {
          bool isProton = isSelectedProton(track);
          if (isProton) {
            float relativeMomentum = computeRelativeMomentum(track, pVec2Prong, massD0);
            if (relativeMomentum < femtoMaxRelativeMomentum) {
              keepEvent[kFemto] = true;
              if (activateQA) {
                hCharmProtonKstarDistr[kD0]->Fill(relativeMomentum);
              }
            }
          }
        } // end femto selection

      } // end loop over tracks
    }   // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      bool isDPlus = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DPlusToPiKPi);
      bool isDs = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DsToPiKK);
      bool isLc = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::LcToPKPi);
      bool isXic = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::XicToPKPi);
      if (!isDPlus && !isDs && !isLc && !isXic) { // check if it's a D+, Ds+ or Lc+
        continue;
      }

      auto trackFirst = cand3Prong.index0_as<BigTracksWithProtonPID>();
      auto trackSecond = cand3Prong.index1_as<BigTracksWithProtonPID>();
      auto trackThird = cand3Prong.index2_as<BigTracksWithProtonPID>();
      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      float sign3Prong = trackFirst.signed1Pt() * trackSecond.signed1Pt() * trackThird.signed1Pt();

      // TODO: add ML selections here
      n3Prongs++;

      auto pVec3Prong = RecoDecay::PVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::Pt(pVec3Prong);

      bool isDPlusInMass = false;
      int isDsInMass = 0;
      int isLcInMass = 0;
      int isXicInMass = 0;
      if (isDPlus) {
        isDPlusInMass = isSelectedDplusInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
      }
      if (isDs) {
        isDsInMass = isSelectedDsInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
      }
      if (isLc) {
        isLcInMass = isSelectedLcInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
      }
      if (isXic) {
        isXicInMass = isSelectedXicInMassRange(pVecFirst, pVecThird, pVecSecond, pt3Prong);
      }

      if (pt3Prong >= pTThreshold3Prong) {
        keepEvent[kHighPt] = true;
        if (activateQA) {
          if (isDPlus) {
            hCharmHighPt[kDplus]->Fill(pt3Prong);
          }
          if (isDs) {
            hCharmHighPt[kDs]->Fill(pt3Prong);
          }
          if (isLc) {
            hCharmHighPt[kLc]->Fill(pt3Prong);
          }
          if (isXic) {
            hCharmHighPt[kXic]->Fill(pt3Prong);
          }
        }
      } // end high-pT selection

      for (const auto& track : tracks) { // start loop over tracks

        if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecFourth = {track.px(), track.py(), track.pz()};

        bool specieCharmHypos[kNBeautyParticles - 2] = {isDPlus, isDs, isLc, isXic};
        int specieCharmHyposInMass[kNBeautyParticles - 2] = {(int)isDPlusInMass, isDsInMass, isLcInMass, isXicInMass};
        float massCharmHypos[kNBeautyParticles - 2] = {massDPlus, massDs, massLc, massXic};
        float massBeautyHypos[kNBeautyParticles - 2] = {massB0, massBs, massLb, massXib};
        float deltaMassHypos[kNBeautyParticles - 2] = {deltaMassB0, deltaMassBs, deltaMassLb, deltaMassXib};
        if (track.signed1Pt() * sign3Prong < 0 && isSelectedTrackForBeauty(track, kBeauty4Prong)) {
          for (int iHypo{0}; iHypo < kNBeautyParticles - 2 && !keepEvent[kBeauty]; ++iHypo) {
            if ((iHypo != 1 && specieCharmHyposInMass[iHypo] > 0) || (iHypo == 1 && ((TESTBIT(specieCharmHyposInMass[iHypo], 0) && TESTBIT(specieCharmHyposInMass[iHypo], 2)) || (TESTBIT(specieCharmHyposInMass[iHypo], 1) && TESTBIT(specieCharmHyposInMass[iHypo], 3))))) {
              auto massCandB = RecoDecay::M(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
              if (std::abs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
                keepEvent[kBeauty] = true;
                if (activateQA) {
                  auto pVecBeauty4Prong = RecoDecay::PVec(pVec3Prong, pVecFourth);
                  auto ptCandBeauty4Prong = RecoDecay::Pt(pVecBeauty4Prong);
                  hMassVsPtB[iHypo + 2]->Fill(ptCandBeauty4Prong, massCandB);
                }
              }
            }
          }
        } // end beauty selection

        // 3-prong femto
        if (isSelectedProton(track)) {
          for (int iHypo{0}; iHypo < kNCharmParticles - 1 && !keepEvent[kFemto]; ++iHypo) {
            if ((iHypo != 1 && specieCharmHypos[iHypo]) || (iHypo == 1 && specieCharmHypos[iHypo] && (TESTBIT(specieCharmHyposInMass[iHypo], 2) || TESTBIT(specieCharmHyposInMass[iHypo], 3)))) {
              float relativeMomentum = computeRelativeMomentum(track, pVec3Prong, massCharmHypos[iHypo]);
              if (relativeMomentum < femtoMaxRelativeMomentum) {
                keepEvent[kFemto] = true;
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

    if (n2Prongs > 1 || n3Prongs > 1 || (n2Prongs > 0 && n3Prongs > 0)) {
      keepEvent[kDoubleCharm] = true;
    }

    tags(keepEvent[kHighPt], keepEvent[kBeauty], keepEvent[kFemto], keepEvent[kDoubleCharm]);

    if (!keepEvent[kHighPt] && !keepEvent[kBeauty] && !keepEvent[kFemto] && !keepEvent[kDoubleCharm]) {
      hProcessedEvents->Fill(1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersHF; ++iTrigger) {
        if (keepEvent[iTrigger]) {
          hProcessedEvents->Fill(iTrigger + 2);
        }
      }
    }
  }

  void processTraining(aod::Hf2Prong const& cand2Prongs,
                       aod::Hf3Prong const& cand3Prongs,
                       aod::McParticles const& particlesMC,
                       BigTracksMCPID const&)
  {
    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      auto trackPos = cand2Prong.index0_as<BigTracksMCPID>(); // positive daughter
      auto trackNeg = cand2Prong.index1_as<BigTracksMCPID>(); // negative daughter

      int8_t sign = 0;
      int8_t flag = 0;
      int8_t origin = 0;

      // D0(bar) → π± K∓
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, std::array{trackPos, trackNeg}, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = (RecoDecay::getMother(particlesMC, particle, kBottom, true) > -1 ? OriginType::NonPrompt : OriginType::Prompt);
        if (origin == OriginType::NonPrompt) {
          flag = kNonPrompt;
        } else {
          flag = kPrompt;
        }
      } else {
        flag = kBkg;
      }

      double pseudoRndm = trackPos.pt() * 1000. - (long)(trackPos.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train2P(trackPos.pt(), trackPos.dcaXY(), trackPos.dcaZ(), trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), trackPos.tofNSigmaPi(), trackPos.tofNSigmaKa(),
                trackNeg.pt(), trackNeg.dcaXY(), trackNeg.dcaZ(), trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaPi(), trackNeg.tofNSigmaKa(),
                flag);
      }
    } // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      auto trackFirst = cand3Prong.index0_as<BigTracksMCPID>();  // first daughter
      auto trackSecond = cand3Prong.index1_as<BigTracksMCPID>(); // second daughter
      auto trackThird = cand3Prong.index2_as<BigTracksMCPID>();  // third daughter
      auto arrayDaughters = std::array{trackFirst, trackSecond, trackThird};

      int8_t sign = 0;
      int8_t flag = 0;
      int8_t channel = -1;
      int8_t origin = 0;

      // D± → π± K∓ π±
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign);
      if (indexRec >= 0) {
        channel = kDplus;
      }
      if (indexRec < 0) {
        // Ds± → K± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 431, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign); // TODO: replace hard coded pdg code
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
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign);
        if (indexRec >= 0) {
          channel = kXic;
        }
      }

      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        origin = (RecoDecay::getMother(particlesMC, particle, kBottom, true) > -1 ? OriginType::NonPrompt : OriginType::Prompt);
        if (origin == OriginType::NonPrompt) {
          flag = kNonPrompt;
        } else {
          flag = kPrompt;
        }
      } else {
        flag = kBkg;
      }

      double pseudoRndm = trackFirst.pt() * 1000. - (long)(trackFirst.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < donwSampleBkgFactor)) {
        train3P(trackFirst.pt(), trackFirst.dcaXY(), trackFirst.dcaZ(), trackFirst.tpcNSigmaPi(), trackFirst.tpcNSigmaKa(), trackFirst.tpcNSigmaPr(), trackFirst.tofNSigmaPi(), trackFirst.tofNSigmaKa(), trackFirst.tofNSigmaPr(),
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
