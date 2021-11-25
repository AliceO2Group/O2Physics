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

static const std::vector<std::string> HfTriggerNames{"highPt", "beauty", "femto", "doubleCharm"};

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
  kXic
};

static const float massPi = RecoDecay::getMassPDG(211);
static const float massK = RecoDecay::getMassPDG(321);
static const float massProton = RecoDecay::getMassPDG(2212);
static const float massD0 = RecoDecay::getMassPDG(421);
static const float massDPlus = RecoDecay::getMassPDG(411);
static const float massDs = RecoDecay::getMassPDG(431);
static const float massLc = RecoDecay::getMassPDG(4122);
static const float massBPlus = RecoDecay::getMassPDG(511);
static const float massB0 = RecoDecay::getMassPDG(521);
static const float massBs = RecoDecay::getMassPDG(531);
static const float massLb = RecoDecay::getMassPDG(5122);

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

  void process(aod::HfTrackIndexProng2 const& cand2Prongs,
               aod::HfTrackIndexProng3 const& cand3Prongs,
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

  // parameters for high-pT triggers
  Configurable<float> pTThreshold2Prong{"pTThreshold2Prong", 5., "pT treshold for high pT 2-prong candidates for kHighPt triggers in GeV/c"};
  Configurable<float> pTThreshold3Prong{"pTThreshold3Prong", 5., "pT treshold for high pT 3-prong candidates for kHighPt triggers in GeV/c"};

  // parameters for beauty triggers
  Configurable<float> deltaMassBPlus{"deltaMassBPlus", 0.3, "invariant-mass delta with respect to the B+ mass"};
  Configurable<float> deltaMassB0{"deltaMassB0", 0.3, "invariant-mass delta with respect to the B0 mass"};
  Configurable<float> deltaMassBs{"deltaMassBs", 0.3, "invariant-mass delta with respect to the Bs mass"};
  Configurable<float> deltaMassLb{"deltaMassLb", 0.3, "invariant-mass delta with respect to the Lb mass"};
  Configurable<float> pTMinBeautyBachelor{"pTMinBeautyBachelor", 0.5, "minumum pT for bachelor pion track used to build b-hadron candidates"};
  Configurable<std::vector<double>> pTBinsTrack{"pTBinsTrack", std::vector<double>{pTBinsTrack_v}, "track pT bin limits for DCAXY pT-depentend cut"};
  Configurable<LabeledArray<double>> cutsTrackBeauty3Prong{"cutsTrackBeauty3Prong", {cutsTrack[0], npTBinsTrack, nCutVarsTrack, pTBinLabelsTrack, cutVarLabelsTrack}, "Single-track selections per pT bin for 3-prong beauty candidates"};
  Configurable<LabeledArray<double>> cutsTrackBeauty4Prong{"cutsTrackBeauty4Prong", {cutsTrack[0], npTBinsTrack, nCutVarsTrack, pTBinLabelsTrack, cutVarLabelsTrack}, "Single-track selections per pT bin for 4-prong beauty candidates"};
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

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  std::shared_ptr<TH1> hProcessedEvents;

  void init(o2::framework::InitContext&)
  {

    cutsSingleTrackBeauty = {cutsTrackBeauty3Prong, cutsTrackBeauty4Prong};

    hProcessedEvents = registry.add<TH1>("fProcessedEvents", "HF - event filtered;;events", HistType::kTH1F, {{5, -0.5, 4.5}});
    std::array<std::string, 5> eventTitles = {"all", "rejected", "w/ high-#it{p}_{T} candidate", "w/ beauty candidate", "w/ femto candidate"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      hProcessedEvents->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }
  }

  /// Single-track cuts for 2-prongs or 3-prongs
  /// \param track is a track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedTrack(const T& track, const int candType)
  {
    auto pT = track.pt();
    if (pT < pTMinBeautyBachelor) {
      return false;
    }
    auto pTBinTrack = findBin(pTBinsTrack, pT);
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(track.dcaXY()) < cutsSingleTrackBeauty[candType].get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; //minimum DCAxy
    }
    if (std::abs(track.dcaXY()) > cutsSingleTrackBeauty[candType].get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; //maximum DCAxy
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

    return true;

  } //bool isSelectedProton(const T& track)

  /// Computation of the relative momentum between particle pairs
  /// \param track is a track
  /// \param ProtonMass is the mass of a proton
  /// \param CharmCandMomentum is the three momentum of a charm candidate
  /// \param CharmMass is the mass of the charm hadron
  /// \return relative momentum of pair
  template <typename T> //template <typename T, typename C>
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
  } //float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

  using HfTrackIndexProng2withColl = soa::Join<aod::HfTrackIndexProng2, aod::Colls2Prong>;
  using HfTrackIndexProng3withColl = soa::Join<aod::HfTrackIndexProng3, aod::Colls3Prong>;
  using BigTracksWithProtonPID = soa::Join<aod::BigTracksExtended, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using BigTracksMCPID = soa::Join<aod::BigTracksPIDExtended, aod::BigTracksMC>;

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

      // TODO: add ML selections here
      n2Prongs++;

      auto pVec2Prong = RecoDecay::PVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::Pt(pVec2Prong);
      if (pt2Prong >= pTThreshold2Prong) {
        keepEvent[kHighPt] = true;
      }

      for (const auto& track : tracks) { // start loop over tracks

        if (track.globalIndex() == trackPos.globalIndex() || track.globalIndex() == trackNeg.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecThird = {track.px(), track.py(), track.pz()};

        if (!keepEvent[kBeauty]) {
          if (isSelectedTrack(track, kBeauty3Prong)) {                                                    // TODO: add more single track cuts
            auto massCandB = RecoDecay::M(std::array{pVec2Prong, pVecThird}, std::array{massD0, massPi}); // TODO: retrieve D0-D0bar hypothesis to pair with proper signed track
            if (std::abs(massCandB - massBPlus) <= deltaMassBPlus) {
              keepEvent[kBeauty] = true;
            }
          }
        }

        // 2-prong femto
        if (!keepEvent[kFemto]) {
          bool isProton = isSelectedProton(track);
          if (isProton) {
            float RelativeMomentum = computeRelativeMomentum(track, pVec2Prong, massD0);
            if (RelativeMomentum < femtoMaxRelativeMomentum) {
              keepEvent[kFemto] = true;
            }
          }
        } //if(!keepEvent[kFemto])

      } // end loop over tracks
    }   // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      bool isDPlus = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DPlusToPiKPi);
      bool isDs = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DsToPiKK);
      bool isLc = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::LcToPKPi);
      if (!isDPlus && !isDs && !isLc) { // check if it's a D+, Ds+ or Lc+
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
      if (pt3Prong >= pTThreshold3Prong) {
        keepEvent[kHighPt] = true;
      }

      for (const auto& track : tracks) { // start loop over tracks

        if (track.globalIndex() == trackFirst.globalIndex() || track.globalIndex() == trackSecond.globalIndex() || track.globalIndex() == trackThird.globalIndex()) {
          continue;
        }

        std::array<float, 3> pVecFourth = {track.px(), track.py(), track.pz()};

        bool specieCharmHypos[3] = {isDPlus, isDs, isLc};
        float massCharmHypos[3] = {massDPlus, massDs, massLc};
        float massBeautyHypos[3] = {massB0, massBs, massLb};
        float deltaMassHypos[3] = {deltaMassB0, deltaMassBs, deltaMassLb};
        if (track.signed1Pt() * sign3Prong < 0 && isSelectedTrack(track, kBeauty4Prong)) { // TODO: add more single track cuts
          for (int iHypo{0}; iHypo < 3 && !keepEvent[kBeauty]; ++iHypo) {
            if (specieCharmHypos[iHypo]) {
              auto massCandB = RecoDecay::M(std::array{pVec3Prong, pVecFourth}, std::array{massCharmHypos[iHypo], massPi});
              if (std::abs(massCandB - massBeautyHypos[iHypo]) <= deltaMassHypos[iHypo]) {
                keepEvent[kBeauty] = true;
              }
            }
          }
        }

        // 3-prong femto
        if (isSelectedProton(track)) {
          for (int iHypo{0}; iHypo < 3 && !keepEvent[kFemto]; ++iHypo) {
            if (specieCharmHypos[iHypo]) {
              float RelativeMomentum = computeRelativeMomentum(track, pVec3Prong, massCharmHypos[iHypo]);
              if (RelativeMomentum < femtoMaxRelativeMomentum) {
                keepEvent[kFemto] = true;
              }
            }
          }
        }
      } // end loop over tracks
    }   // end loop over 3-prong candidates

    if (n2Prongs > 1 || n3Prongs > 1 || (n2Prongs > 0 && n3Prongs > 0)) {
      keepEvent[kDoubleCharm] = true;
    }

    tags(keepEvent[kHighPt], keepEvent[kBeauty], keepEvent[kFemto], keepEvent[kDoubleCharm]);

    if (!keepEvent[kHighPt] && !keepEvent[kBeauty] && !keepEvent[kFemto] && keepEvent[kDoubleCharm]) {
      hProcessedEvents->Fill(1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersHF; iTrigger++) {
        if (keepEvent[iTrigger]) {
          hProcessedEvents->Fill(iTrigger + 2);
        }
      }
    }
  }

  void processTraining(aod::HfTrackIndexProng2 const& cand2Prongs,
                       aod::HfTrackIndexProng3 const& cand3Prongs,
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
        auto particle = particlesMC.iteratorAt(indexRec);
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
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 431, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign); //TODO: replace hard coded pdg code
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
        auto particle = particlesMC.iteratorAt(indexRec);
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
