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

/// \author Jasper Parkkila <jasper.parkkila@cern.ch>

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/DataModel/PIDResponseITS.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/detail/TypeTruncation.h"

#include <TFormula.h>

#include <experimental/type_traits>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

enum LambdaPid { kLambda = 0,
                 kAntiLambda
};

// #define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct Filter2Prong {
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 0, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgYMax, float, -1.0f, "Maximum candidate rapidity")
  O2_DEFINE_CONFIGURABLE(cfgImPart1Mass, float, o2::constants::physics::MassKPlus, "Daughter particle 1 mass in GeV")
  O2_DEFINE_CONFIGURABLE(cfgImPart2Mass, float, o2::constants::physics::MassKMinus, "Daughter particle 2 mass in GeV")
  O2_DEFINE_CONFIGURABLE(cfgImPart1PID, int, o2::track::PID::Kaon, "PID of daughter particle 1 (O2 PID ID)")
  O2_DEFINE_CONFIGURABLE(cfgImPart2PID, int, o2::track::PID::Kaon, "PID of daughter particle 2 (O2 PID ID)")
  O2_DEFINE_CONFIGURABLE(cfgMomDepPID, bool, 1, "Use mommentum dependent PID for Phi meson")
  O2_DEFINE_CONFIGURABLE(cfgImCutPt, float, 0.2f, "Minimal pT for candidates")
  O2_DEFINE_CONFIGURABLE(cfgImMinInvMass, float, 0.95f, "Minimum invariant mass for generic 2 prong")
  O2_DEFINE_CONFIGURABLE(cfgImMaxInvMass, float, 1.07f, "Maximum invariant mass for generic 2 prong")
  O2_DEFINE_CONFIGURABLE(cfgImSigmaFormula, std::string, "(([p] < 0.5 || [hasTOF] <= 0.0) && abs([sTPC]) < 3.0) || ([p] >= 0.5 && abs([sTPC]) < 2.5 && abs([sTOF]) < 3.0)", "pT dependent daughter track sigma pass condition. Parameters: [p] momentum, [sTPC] sigma TPC, [sTOF] sigma TOF, [hasTOF] has TOF.")

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(tpcNClsCrossedRowsTrackMin, float, 70, "Minimum number of crossed rows in TPC");
    O2_DEFINE_CONFIGURABLE(etaTrackMax, float, 0.8, "Maximum pseudorapidity");
    O2_DEFINE_CONFIGURABLE(ptTrackMin, float, 0.1, "Minimum transverse momentum");
    O2_DEFINE_CONFIGURABLE(minV0DCAPr, float, 0.1, "Min V0 proton DCA");
    O2_DEFINE_CONFIGURABLE(minV0DCAPiLambda, float, 0.1, "Min V0 pion DCA for lambda");
    O2_DEFINE_CONFIGURABLE(minV0DCAPiK0s, float, 0.1, "Min V0 pion DCA for K0s");
    O2_DEFINE_CONFIGURABLE(daughPIDCuts, float, 4.0, "PID nsigma for V0s");
    O2_DEFINE_CONFIGURABLE(massK0Min, float, 0.4, "Minimum mass for K0");
    O2_DEFINE_CONFIGURABLE(massK0Max, float, 0.6, "Maximum mass for K0");
    O2_DEFINE_CONFIGURABLE(massLambdaMin, float, 1.0, "Minimum mass for lambda");
    O2_DEFINE_CONFIGURABLE(massLambdaMax, float, 1.3, "Maximum mass for lambda");
    O2_DEFINE_CONFIGURABLE(radiusMaxLambda, float, 2.3, "Maximum decay radius (cm) for lambda");
    O2_DEFINE_CONFIGURABLE(radiusMinLambda, float, 0.0, "Minimum decay radius (cm) for lambda");
    O2_DEFINE_CONFIGURABLE(radiusMaxK0s, float, 2.3, "Maximum decay radius (cm) for K0s");
    O2_DEFINE_CONFIGURABLE(radiusMinK0s, float, 0.0, "Minimum decay radius (cm) for K0s");
    O2_DEFINE_CONFIGURABLE(cosPaMinLambda, float, 0.98, "Minimum cosine of pointing angle for lambda");
    O2_DEFINE_CONFIGURABLE(cosPaMinK0s, float, 0.98, "Minimum cosine of pointing angle for K0s");
    O2_DEFINE_CONFIGURABLE(dcaV0DaughtersMaxLambda, float, 0.2, "Maximum DCA among the V0 daughters (cm) for lambda");
    O2_DEFINE_CONFIGURABLE(dcaV0DaughtersMaxK0s, float, 0.2, "Maximum DCA among the V0 daughters (cm) for K0s");
    O2_DEFINE_CONFIGURABLE(qtArmenterosMinForK0s, float, 0.12, "Minimum Armenteros' qt for K0s");
    O2_DEFINE_CONFIGURABLE(maxLambdaLifeTime, float, 30, "Maximum lambda lifetime (in cm)");
    O2_DEFINE_CONFIGURABLE(maxK0sLifeTime, float, 30, "Maximum K0s lifetime (in cm)");
  } grpV0;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(ImMinInvMassPhiMeson, float, 0.98f, "Minimum invariant mass Phi meson (GeV)");
    O2_DEFINE_CONFIGURABLE(ImMaxInvMassPhiMeson, float, 1.07f, "Maximum invariant mass Phi meson (GeV)");
    O2_DEFINE_CONFIGURABLE(ITSPIDSelection, bool, true, "PID ITS");
    O2_DEFINE_CONFIGURABLE(ITSPIDPthreshold, float, 1.0, "Momentum threshold for ITS PID (GeV/c) (only used if ITSPIDSelection is true)");
    O2_DEFINE_CONFIGURABLE(lowITSPIDNsigma, float, 3.0, "lower cut on PID nsigma for ITS");
    O2_DEFINE_CONFIGURABLE(highITSPIDNsigma, float, 3.0, "higher cut on PID nsigma for ITS");
    O2_DEFINE_CONFIGURABLE(ITSclusterPhiMeson, int, 5, "Minimum number of ITS cluster for phi meson track");
    O2_DEFINE_CONFIGURABLE(TPCCrossedRowsPhiMeson, int, 80, "Minimum number of TPC Crossed Rows for phi meson track");
    O2_DEFINE_CONFIGURABLE(cutDCAxyPhiMeson, float, 0.1, "Maximum DCAxy for phi meson track");
    O2_DEFINE_CONFIGURABLE(cutDCAzPhiMeson, float, 0.1, "Maximum DCAz for phi meson track");
    O2_DEFINE_CONFIGURABLE(cutEtaPhiMeson, float, 0.8, "Maximum eta for phi meson track");
    O2_DEFINE_CONFIGURABLE(cutPTPhiMeson, float, 0.15, "Maximum pt for phi meson track");
    O2_DEFINE_CONFIGURABLE(isDeepAngle, bool, true, "Flag for applying deep angle");
    O2_DEFINE_CONFIGURABLE(deepAngle, float, 0.04, "Deep angle cut");
    O2_DEFINE_CONFIGURABLE(nsigmaCutTPC, float, 2.5, "nsigma tpc");
    O2_DEFINE_CONFIGURABLE(nsigmaCutTOF, float, 2.5, "nsigma tof");
    O2_DEFINE_CONFIGURABLE(cutTOFBeta, float, 0.5, "TOF beta");
    O2_DEFINE_CONFIGURABLE(confFakeKaonCut, float, 0.15, "Cut based on track from momentum difference");
    O2_DEFINE_CONFIGURABLE(removefaketrack, bool, true, "Flag to remove fake kaon");
  } grpPhi;

  HfHelper hfHelper;
  Produces<aod::CF2ProngTracks> output2ProngTracks;
  Produces<aod::CF2ProngTrackmls> output2ProngTrackmls;

  Produces<aod::CF2ProngMcParts> output2ProngMcParts;

  std::vector<float> mlvecd{};
  std::vector<float> mlvecdbar{};

  using HFCandidates = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using HFCandidatesML = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;

  template <class T>
  using HasMLProb = decltype(std::declval<T&>().mlProbD0());

  using PIDTrack = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta, aod::TracksDCA>;
  using ResoV0s = aod::V0Datas;

  std::unique_ptr<TFormula> sigmaFormula;
  std::array<uint, 4> sigmaFormulaParamIndex;

  void init(InitContext&)
  {
    if (doprocessDataInvMass) {
      sigmaFormula = std::make_unique<TFormula>("sigmaFormula", cfgImSigmaFormula.value.c_str());
      if (static_cast<std::size_t>(sigmaFormula->GetNpar()) > std::size(sigmaFormulaParamIndex))
        LOGF(fatal, "Number of parameters in cfgImSigmaFormula can not be larger than %d.", std::size(sigmaFormulaParamIndex));
      // could do SetParameter(name,value) directly, but pre-lookup of the names will result in faster process
      std::array<std::string, 4> pars = {"p", "sTPC", "sTOF", "hasTOF"};
      std::fill_n(sigmaFormulaParamIndex.begin(), std::size(sigmaFormulaParamIndex), ~0u);
      for (uint i = 0, n = sigmaFormula->GetNpar(); i < n; ++i) {
        auto m = std::find(pars.begin(), pars.end(), sigmaFormula->GetParName(i));
        if (m != pars.end())
          sigmaFormulaParamIndex[std::distance(pars.begin(), m)] = i;
        else
          LOGF(warning, "Unrecognized cfgImSigmaFormula parameter %s.", sigmaFormula->GetParName(i));
      }
    }
  }

  template <class HFCandidatesType>
  void processDataT(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesType const& candidates)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision
    if (cfgVerbosity > 0 && candidates.size() > 0)
      LOGF(info, "Candidates for collision: %lu, cfcollisions: %lu, CFTracks: %lu", candidates.size(), cfcollisions.size(), cftracks.size());
    for (const auto& c : candidates) {
      int prongCFId[2] = {-1, -1};
      for (const auto& cftrack : cftracks) {
        if (c.prong0Id() == cftrack.trackId()) {
          prongCFId[0] = cftrack.globalIndex();
          break;
        }
      }
      for (const auto& cftrack : cftracks) {
        if (c.prong1Id() == cftrack.trackId()) {
          prongCFId[1] = cftrack.globalIndex();
          break;
        }
      }
      // look-up the collision id
      if ((c.hfflag() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) == 0)
        continue;
      if (cfgYMax >= 0.0f && std::abs(hfHelper.yD0(c)) > cfgYMax)
        continue;

      if (c.isSelD0() > 0) {
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           prongCFId[0], prongCFId[1], c.pt(), c.eta(), c.phi(), hfHelper.invMassD0ToPiK(c), aod::cf2prongtrack::D0ToPiK);
        if constexpr (std::experimental::is_detected<HasMLProb, typename HFCandidatesType::iterator>::value) {
          mlvecd.clear();
          for (const float val : c.mlProbD0())
            mlvecd.push_back(val);
          mlvecdbar.clear();
          for (const float val : c.mlProbD0bar())
            mlvecdbar.push_back(val);
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }

      if (c.isSelD0bar() > 0) {
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           prongCFId[0], prongCFId[1], c.pt(), c.eta(), c.phi(), hfHelper.invMassD0barToKPi(c), aod::cf2prongtrack::D0barToKPi);
        if constexpr (std::experimental::is_detected<HasMLProb, typename HFCandidatesType::iterator>::value) {
          mlvecd.clear();
          for (const float val : c.mlProbD0())
            mlvecd.push_back(val);
          mlvecdbar.clear();
          for (const float val : c.mlProbD0bar())
            mlvecdbar.push_back(val);
          output2ProngTrackmls(cfcollisions.begin().globalIndex(), mlvecd, mlvecdbar);
        }
      }
    }
  }

  void processDataML(aod::Collisions::iterator const& col, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidatesML const& candidates)
  {
    processDataT(col, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processDataML, "Process data D0 candidates with ML", false);

  void processData(aod::Collisions::iterator const& col, aod::BCsWithTimestamps const& bcs, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, HFCandidates const& candidates)
  {
    processDataT(col, bcs, cfcollisions, cftracks, candidates);
  }
  PROCESS_SWITCH(Filter2Prong, processData, "Process data D0 candidates", true);

  using HFMCTrack = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  void processMC(aod::McCollisions::iterator const&, aod::CFMcParticleRefs const& cfmcparticles, [[maybe_unused]] HFMCTrack const& mcparticles)
  {
    // The main filter outputs the primary MC particles. Here we just resolve the daughter indices that are needed for the efficiency matching.
    for (const auto& r : cfmcparticles) {
      const auto& mcParticle = r.mcParticle_as<HFMCTrack>();
      if ((std::abs(mcParticle.flagMcMatchGen()) != o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) || mcParticle.daughtersIds().size() != 2) {
        output2ProngMcParts(-1, -1, aod::cf2prongtrack::Generic2Prong);
        continue;
      }
      int prongCFId[2] = {-1, -1};
      for (uint i = 0; i < 2; ++i) {
        for (const auto& cfmcpart : cfmcparticles) {
          if (mcParticle.daughtersIds()[i] == cfmcpart.mcParticleId()) {
            prongCFId[i] = cfmcpart.globalIndex();
            break;
          }
        }
      }
      output2ProngMcParts(prongCFId[0], prongCFId[1],
                          (mcParticle.pdgCode() >= 0 ? aod::cf2prongtrack::D0ToPiK : aod::cf2prongtrack::D0barToKPi) | ((mcParticle.originMcGen() & RecoDecay::OriginType::Prompt) ? aod::cf2prongmcpart::Prompt : 0));
    }
  }
  PROCESS_SWITCH(Filter2Prong, processMC, "Process MC 2-prong daughters", false);

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() >= grpPhi.ITSclusterPhiMeson && candidate.tpcNClsCrossedRows() > grpPhi.TPCCrossedRowsPhiMeson && std::abs(candidate.dcaXY()) <= grpPhi.cutDCAxyPhiMeson && std::abs(candidate.dcaZ()) <= grpPhi.cutDCAzPhiMeson && std::abs(candidate.eta()) <= grpPhi.cutEtaPhiMeson && candidate.pt() >= grpPhi.cutPTPhiMeson) {
      return true;
    }
    return false;
  }

  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (grpPhi.isDeepAngle && angle < grpPhi.deepAngle) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isFakeTrack(T const& track)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (TMath::Abs(pglobal - ptpc) > grpPhi.confFakeKaonCut) {
      return true;
    }
    return false;
  }

  template <typename Collision, typename V0Cand>
  bool isSelectedV0AsK0s(Collision const& collision, const V0Cand& v0)
  {
    const auto& posTrack = v0.template posTrack_as<PIDTrack>();
    const auto& negTrack = v0.template negTrack_as<PIDTrack>();

    float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0;

    if (v0.mK0Short() < grpV0.massK0Min || v0.mK0Short() > grpV0.massK0Max) {
      return false;
    }
    if ((v0.qtarm() / std::abs(v0.alpha())) < grpV0.qtArmenterosMinForK0s) {
      return false;
    }
    if (v0.v0radius() > grpV0.radiusMaxK0s || v0.v0radius() < grpV0.radiusMinK0s) {
      return false;
    }
    if (v0.v0cosPA() < grpV0.cosPaMinK0s) {
      return false;
    }
    if (v0.dcaV0daughters() > grpV0.dcaV0DaughtersMaxK0s) {
      return false;
    }
    if (std::abs(CtauK0s) > grpV0.maxK0sLifeTime) {
      return false;
    }
    if (((std::abs(posTrack.tpcNSigmaPi()) > grpV0.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > grpV0.daughPIDCuts))) {
      return false;
    }
    if ((TMath::Abs(v0.dcapostopv()) < grpV0.minV0DCAPiK0s || TMath::Abs(v0.dcanegtopv()) < grpV0.minV0DCAPiK0s)) {
      return false;
    }
    return true;
  }

  template <LambdaPid pid, typename Collision, typename V0Cand>
  bool isSelectedV0AsLambda(Collision const& collision, const V0Cand& v0)
  {
    const auto& posTrack = v0.template posTrack_as<PIDTrack>();
    const auto& negTrack = v0.template negTrack_as<PIDTrack>();

    float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda;

    if ((v0.mLambda() < grpV0.massLambdaMin || v0.mLambda() > grpV0.massLambdaMax) &&
        (v0.mAntiLambda() < grpV0.massLambdaMin || v0.mAntiLambda() > grpV0.massLambdaMax)) {
      return false;
    }
    if (v0.v0radius() > grpV0.radiusMaxLambda || v0.v0radius() < grpV0.radiusMinLambda) {
      return false;
    }
    if (v0.v0cosPA() < grpV0.cosPaMinLambda) {
      return false;
    }
    if (v0.dcaV0daughters() > grpV0.dcaV0DaughtersMaxLambda) {
      return false;
    }
    if (pid == LambdaPid::kLambda && (TMath::Abs(v0.dcapostopv()) < grpV0.minV0DCAPr || TMath::Abs(v0.dcanegtopv()) < grpV0.minV0DCAPiLambda)) {
      return false;
    }
    if (pid == LambdaPid::kAntiLambda && (TMath::Abs(v0.dcapostopv()) < grpV0.minV0DCAPiLambda || TMath::Abs(v0.dcanegtopv()) < grpV0.minV0DCAPr)) {
      return false;
    }
    if (pid == LambdaPid::kLambda && ((std::abs(posTrack.tpcNSigmaPr()) > grpV0.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > grpV0.daughPIDCuts))) {
      return false;
    }
    if (pid == LambdaPid::kAntiLambda && ((std::abs(posTrack.tpcNSigmaPi()) > grpV0.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPr()) > grpV0.daughPIDCuts))) {
      return false;
    }
    if (std::abs(CtauLambda) > grpV0.maxLambdaLifeTime) {
      return false;
    }
    return true;
  }

  template <typename T1>
  bool isV0TrackSelected(const T1& v0)
  {
    const auto& posTrack = v0.template posTrack_as<PIDTrack>();
    const auto& negTrack = v0.template negTrack_as<PIDTrack>();

    if (!posTrack.hasTPC() || !negTrack.hasTPC()) {
      return false;
    }
    if (posTrack.tpcNClsCrossedRows() < grpV0.tpcNClsCrossedRowsTrackMin || negTrack.tpcNClsCrossedRows() < grpV0.tpcNClsCrossedRowsTrackMin) {
      return false;
    }
    if (posTrack.tpcCrossedRowsOverFindableCls() < 0.8 || negTrack.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }
    if (std::abs(v0.positiveeta()) > grpV0.etaTrackMax || std::abs(v0.negativeeta()) > grpV0.etaTrackMax) {
      return false;
    }
    if (v0.positivept() < grpV0.ptTrackMin || v0.negativept() < grpV0.ptTrackMin) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (cfgMomDepPID) {
      if (candidate.p() < 0.5) {
        if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < grpPhi.nsigmaCutTPC) {
          return true;
        } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < grpPhi.nsigmaCutTOF && candidate.beta() > grpPhi.cutTOFBeta) {
          return true;
        }
      } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < grpPhi.nsigmaCutTOF && candidate.beta() > grpPhi.cutTOFBeta) {
        return true;
      }
    } else {
      if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < grpPhi.nsigmaCutTPC) {
        return true;
      } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < grpPhi.nsigmaCutTOF && candidate.beta() > grpPhi.cutTOFBeta) {
        return true;
      }
    }
    return false;
  }

  // Generic 2-prong invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  void processDataInvMass(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, Filter2Prong::PIDTrack const& tracks)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision
    for (const auto& cftrack1 : cftracks) {
      const auto& p1 = tracks.iteratorAt(cftrack1.trackId() - tracks.begin().globalIndex());
      if (p1.sign() != 1)
        continue;
      sigmaFormula->SetParameter(sigmaFormulaParamIndex[0], p1.p());
      sigmaFormula->SetParameter(sigmaFormulaParamIndex[1], o2::aod::pidutils::tpcNSigma(cfgImPart1PID, p1));
      sigmaFormula->SetParameter(sigmaFormulaParamIndex[2], o2::aod::pidutils::tofNSigma(cfgImPart1PID, p1));
      sigmaFormula->SetParameter(sigmaFormulaParamIndex[3], static_cast<double>(p1.hasTOF()));
      if (sigmaFormula->Eval() <= 0.0f)
        continue;
      for (const auto& cftrack2 : cftracks) {
        if (cftrack2.globalIndex() == cftrack1.globalIndex())
          continue;
        const auto& p2 = tracks.iteratorAt(cftrack2.trackId() - tracks.begin().globalIndex());
        if (p2.sign() != -1)
          continue;
        sigmaFormula->SetParameter(sigmaFormulaParamIndex[0], p2.p());
        sigmaFormula->SetParameter(sigmaFormulaParamIndex[1], o2::aod::pidutils::tpcNSigma(cfgImPart1PID, p2));
        sigmaFormula->SetParameter(sigmaFormulaParamIndex[2], o2::aod::pidutils::tofNSigma(cfgImPart1PID, p2));
        sigmaFormula->SetParameter(sigmaFormulaParamIndex[3], static_cast<double>(p2.hasTOF()));
        if (sigmaFormula->Eval() <= 0.0f)
          continue;
        ROOT::Math::PtEtaPhiMVector vec1(p1.pt(), p1.eta(), p1.phi(), cfgImPart1Mass);
        ROOT::Math::PtEtaPhiMVector vec2(p2.pt(), p2.eta(), p2.phi(), cfgImPart2Mass);
        ROOT::Math::PtEtaPhiMVector s = vec1 + vec2;
        if (s.pt() < cfgImCutPt || s.M() < cfgImMinInvMass || s.M() > cfgImMaxInvMass)
          continue;

        float phi = RecoDecay::constrainAngle(s.Phi(), 0.0f);
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           cftrack1.globalIndex(), cftrack2.globalIndex(), s.pt(), s.eta(), phi, s.M(), aod::cf2prongtrack::Generic2Prong);
      }
    }
  }
  PROCESS_SWITCH(Filter2Prong, processDataInvMass, "Process data generic 2-prong candidates with invariant mass method", false);

  // Phi and V0s invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  void processDataPhiV0(aod::Collisions::iterator const& collision, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, Filter2Prong::PIDTrack const& tracks, aod::V0Datas const& V0s)
  {
    if (cfcollisions.size() <= 0)
      return; // rejected collision

    // V0
    for (const auto& v0 : V0s) {    // Loop over V0 candidates
      if (!isV0TrackSelected(v0)) { // Quality selection for V0 prongs
        continue;
      }

      const auto& posTrack = v0.template posTrack_as<PIDTrack>();
      const auto& negTrack = v0.template negTrack_as<PIDTrack>();
      double massV0 = 0.0;

      // K0s
      if (isSelectedV0AsK0s(collision, v0)) { // candidate is K0s
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), v0.mK0Short(), aod::cf2prongtrack::K0stoPiPi);
      }

      // Lambda and Anti-Lambda
      bool LambdaTag = isSelectedV0AsLambda<LambdaPid::kLambda>(collision, v0);
      bool aLambdaTag = isSelectedV0AsLambda<LambdaPid::kAntiLambda>(collision, v0);

      // Note: candidate compatible with Lambda and Anti-Lambda hypothesis are counted twice (once for each hypothesis)
      if (LambdaTag) { // candidate is Lambda
        massV0 = v0.mLambda();
        output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), massV0, aod::cf2prongtrack::LambdatoPPi);
      }
      if (aLambdaTag) { // candidate is Anti-lambda
        massV0 = v0.mAntiLambda();
        output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), massV0, aod::cf2prongtrack::AntiLambdatoPiP);
      } // end of Lambda and Anti-Lambda processing
    } // end of loop over V0 candidates

    // Phi
    if (cftracks.size() <= 0)
      return; // rejected collision

    o2::aod::ITSResponse itsResponse;

    for (const auto& cftrack1 : cftracks) { // Loop over first track
      const auto& p1 = tracks.iteratorAt(cftrack1.trackId() - tracks.begin().globalIndex());
      if (p1.sign() != 1) {
        continue;
      }
      if (!selectionTrack(p1)) {
        continue;
      }
      if (grpPhi.ITSPIDSelection && p1.p() < grpPhi.ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) > grpPhi.lowITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) < grpPhi.highITSPIDNsigma.value)) { // Check ITS PID condition
        continue;
      }
      if (!selectionPID(p1)) {
        continue;
      }
      if (grpPhi.removefaketrack && isFakeTrack(p1)) { // Check if the track is a fake kaon
        continue;
      }

      for (const auto& cftrack2 : cftracks) {                 // Loop over second track
        if (cftrack2.globalIndex() == cftrack1.globalIndex()) // Skip if it's the same track as the first one
          continue;

        const auto& p2 = tracks.iteratorAt(cftrack2.trackId() - tracks.begin().globalIndex());
        if (p2.sign() != -1) {
          continue;
        }
        if (!selectionTrack(p2)) {
          continue;
        }
        if (!selectionPID(p2)) {
          continue;
        }
        if (grpPhi.ITSPIDSelection && p2.p() < grpPhi.ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) > grpPhi.lowITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) < grpPhi.highITSPIDNsigma.value)) { // Check ITS PID condition
          continue;
        }
        if (grpPhi.removefaketrack && isFakeTrack(p2)) { // Check if the track is a fake kaon
          continue;
        }
        if (!selectionPair(p1, p2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector vec1(p1.pt(), p1.eta(), p1.phi(), cfgImPart1Mass);
        ROOT::Math::PtEtaPhiMVector vec2(p2.pt(), p2.eta(), p2.phi(), cfgImPart2Mass);
        ROOT::Math::PtEtaPhiMVector s = vec1 + vec2;
        if (s.M() < grpPhi.ImMinInvMassPhiMeson || s.M() > grpPhi.ImMaxInvMassPhiMeson) {
          continue;
        }
        float phi = RecoDecay::constrainAngle(s.Phi(), 0.0f);
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           cftrack1.globalIndex(), cftrack2.globalIndex(), s.pt(), s.eta(), phi, s.M(), aod::cf2prongtrack::PhiToKK);
      } // end of loop over second track
    } // end of loop over first track
  }
  PROCESS_SWITCH(Filter2Prong, processDataPhiV0, "Process data Phi and V0 candidates with invariant mass method", false);

  // Phi and V0s invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  void processDataV0(aod::Collisions::iterator const& collision, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, Filter2Prong::PIDTrack const&, aod::V0Datas const& V0s)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision

    for (const auto& v0 : V0s) {    // Loop over V0 candidates
      if (!isV0TrackSelected(v0)) { // Quality selection for V0 prongs
        continue;
      }

      const auto& posTrack = v0.template posTrack_as<PIDTrack>();
      const auto& negTrack = v0.template negTrack_as<PIDTrack>();
      double massV0 = 0.0;

      // K0s
      if (isSelectedV0AsK0s(collision, v0)) { // candidate is K0s
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), v0.mK0Short(), aod::cf2prongtrack::K0stoPiPi);
      }

      // Lambda and Anti-Lambda
      bool LambdaTag = isSelectedV0AsLambda<LambdaPid::kLambda>(collision, v0);
      bool aLambdaTag = isSelectedV0AsLambda<LambdaPid::kAntiLambda>(collision, v0);

      // Note: candidate compatible with Lambda and Anti-Lambda hypothesis are counted twice (once for each hypothesis)
      if (LambdaTag) { // candidate is Lambda
        massV0 = v0.mLambda();
        output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), massV0, aod::cf2prongtrack::LambdatoPPi);
      }
      if (aLambdaTag) { // candidate is Anti-lambda
        massV0 = v0.mAntiLambda();
        output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                           v0.pt(), v0.eta(), v0.phi(), massV0, aod::cf2prongtrack::AntiLambdatoPiP);
      } // end of Lambda and Anti-Lambda processing
    } // end of loop over V0 candidates
  }
  PROCESS_SWITCH(Filter2Prong, processDataV0, "Process data V0 candidates with invariant mass method", false);

  // Phi and V0s invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  void processDataPhi(aod::Collisions::iterator const&, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, Filter2Prong::PIDTrack const& tracks)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision

    o2::aod::ITSResponse itsResponse;

    for (const auto& cftrack1 : cftracks) { // Loop over first track
      const auto& p1 = tracks.iteratorAt(cftrack1.trackId() - tracks.begin().globalIndex());
      if (p1.sign() != 1) {
        continue;
      }
      if (!selectionTrack(p1)) {
        continue;
      }
      if (grpPhi.ITSPIDSelection && p1.p() < grpPhi.ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) > grpPhi.lowITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) < grpPhi.highITSPIDNsigma.value)) { // Check ITS PID condition
        continue;
      }
      if (!selectionPID(p1)) {
        continue;
      }
      if (grpPhi.removefaketrack && isFakeTrack(p1)) { // Check if the track is a fake kaon
        continue;
      }

      for (const auto& cftrack2 : cftracks) {                 // Loop over second track
        if (cftrack2.globalIndex() == cftrack1.globalIndex()) // Skip if it's the same track as the first one
          continue;

        const auto& p2 = tracks.iteratorAt(cftrack2.trackId() - tracks.begin().globalIndex());
        if (p2.sign() != -1) {
          continue;
        }
        if (!selectionTrack(p2)) {
          continue;
        }
        if (!selectionPID(p2)) {
          continue;
        }
        if (grpPhi.ITSPIDSelection && p2.p() < grpPhi.ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) > grpPhi.lowITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) < grpPhi.highITSPIDNsigma.value)) { // Check ITS PID condition
          continue;
        }
        if (grpPhi.removefaketrack && isFakeTrack(p2)) { // Check if the track is a fake kaon
          continue;
        }
        if (!selectionPair(p1, p2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector vec1(p1.pt(), p1.eta(), p1.phi(), cfgImPart1Mass);
        ROOT::Math::PtEtaPhiMVector vec2(p2.pt(), p2.eta(), p2.phi(), cfgImPart2Mass);
        ROOT::Math::PtEtaPhiMVector s = vec1 + vec2;
        if (s.M() < grpPhi.ImMinInvMassPhiMeson || s.M() > grpPhi.ImMaxInvMassPhiMeson) {
          continue;
        }
        float phi = RecoDecay::constrainAngle(s.Phi(), 0.0f);
        output2ProngTracks(cfcollisions.begin().globalIndex(),
                           cftrack1.globalIndex(), cftrack2.globalIndex(), s.pt(), s.eta(), phi, s.M(), aod::cf2prongtrack::PhiToKK);
      } // end of loop over second track
    } // end of loop over first track
  }
  PROCESS_SWITCH(Filter2Prong, processDataPhi, "Process data Phi candidates with invariant mass method", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Filter2Prong>(cfgc, TaskName{"filter-2prong"})};
}
