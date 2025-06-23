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

// #define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct Filter2Prong {
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 0, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgYMax, float, -1.0f, "Maximum candidate rapidity")
  O2_DEFINE_CONFIGURABLE(cfgImPart1Mass, float, o2::constants::physics::MassKPlus, "Daughter particle 1 mass in GeV")
  O2_DEFINE_CONFIGURABLE(cfgImPart2Mass, float, o2::constants::physics::MassKMinus, "Daughter particle 2 mass in GeV")
  // O2_DEFINE_CONFIGURABLE(cfgImPart1PID, int, o2::track::PID::Kaon, "PID of daughter particle 1 (O2 PID ID)")
  // O2_DEFINE_CONFIGURABLE(cfgImPart2PID, int, o2::track::PID::Kaon, "PID of daughter particle 2 (O2 PID ID)")
  O2_DEFINE_CONFIGURABLE(cfgMomDepPID, bool, 1, "Use mommentum dependent PID for Phi meson")
  O2_DEFINE_CONFIGURABLE(cfgImCutPt, float, 0.2f, "Minimal pT for candidates")
  O2_DEFINE_CONFIGURABLE(cfgImMinInvMass, float, 0.95f, "Minimum invariant mass (GeV)")
  O2_DEFINE_CONFIGURABLE(cfgImMaxInvMass, float, 1.07f, "Maximum invariant mass (GeV)")
  // O2_DEFINE_CONFIGURABLE(cfgImSigmaFormula, std::string, "(z < 0.5 && x < 3.0) || (z >= 0.5 && x < 2.5 && y < 3.0)", "pT dependent daughter track sigma pass condition (x = TPC sigma, y = TOF sigma, z = pT)")
  O2_DEFINE_CONFIGURABLE(cfgDoPhi, bool, false, "Store phi information")
  O2_DEFINE_CONFIGURABLE(cfgDoV0, bool, true, "Store V0s candidates")
  O2_DEFINE_CONFIGURABLE(tpcNClsCrossedRowsTrackMin, float, 70, "Minimum number of crossed rows in TPC")
  O2_DEFINE_CONFIGURABLE(etaTrackMax, float, 0.8, "Maximum pseudorapidity")
  O2_DEFINE_CONFIGURABLE(ptTrackMin, float, 0.1, "Minimum transverse momentum")
  O2_DEFINE_CONFIGURABLE(cMinV0DCAPr, float, 0.1, "Min V0 proton DCA")
  O2_DEFINE_CONFIGURABLE(cMinV0DCAPi, float, 0.1, "Min V0 pion DCA")
  O2_DEFINE_CONFIGURABLE(ITSPIDSelection, bool, true, "PID ITS")
  O2_DEFINE_CONFIGURABLE(ITSPIDPthreshold, float, 1.0, "Momentum threshold for ITS PID (GeV/c) (only used if ITSPIDSelection is true)")
  O2_DEFINE_CONFIGURABLE(ITSPIDNsigma, float, 3.0, "PID nsigma for ITS")
  O2_DEFINE_CONFIGURABLE(ConfDaughPIDCuts, float, 4.0, "PID nsigma for V0s")
  O2_DEFINE_CONFIGURABLE(massK0Min, float, 0.4, "Minimum mass for K0")
  O2_DEFINE_CONFIGURABLE(massK0Max, float, 0.6, "Maximum mass for K0")
  O2_DEFINE_CONFIGURABLE(massLambdaMin, float, 1.0, "Minimum mass for lambda")
  O2_DEFINE_CONFIGURABLE(massLambdaMax, float, 1.3, "Maximum mass for lambda")
  O2_DEFINE_CONFIGURABLE(massOmegaMin, float, 1.5, "Minimum mass for omega")
  O2_DEFINE_CONFIGURABLE(massOmegaMax, float, 1.8, "Maximum mass for omega")
  O2_DEFINE_CONFIGURABLE(interactionRateMin, float, -1, "Minimum interaction rate (kHz)")
  O2_DEFINE_CONFIGURABLE(interactionRateMax, float, 1.e20, "Maximum interaction rate (kHz)")
  O2_DEFINE_CONFIGURABLE(radiusMax, float, 2.3, "Maximum decay radius (cm)")
  O2_DEFINE_CONFIGURABLE(radiusMin, float, 0.0, "Minimum decay radius (cm)")
  O2_DEFINE_CONFIGURABLE(cosPaMin, float, 0.98, "Minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(dcaV0DaughtersMax, float, 0.2, "Maximum DCA among the V0 daughters (cm)")
  O2_DEFINE_CONFIGURABLE(dcaV0ToPvMax, float, 0.2, "Maximum DCA of the V0 from the primary vertex (cm)")
  O2_DEFINE_CONFIGURABLE(cosPaV0Min, float, 0.95, "Minimum cosine of pointing angle for V0 stemming from cascade decays")
  O2_DEFINE_CONFIGURABLE(qtArmenterosMinForK0, float, 0.12, "Minimum Armenteros' qt for K0")
  O2_DEFINE_CONFIGURABLE(qtArmenterosMaxForLambda, float, 0.12, "Minimum Armenteros' qt for (anti)Lambda")
  O2_DEFINE_CONFIGURABLE(ConfV0Rap, float, 0.5, "Store rapidity of v0")
  O2_DEFINE_CONFIGURABLE(cMaxLambdaLifeTime, float, 30, "Store Lambda lifetime")
  O2_DEFINE_CONFIGURABLE(cMaxK0sLifeTime, float, 30, "Store K0s lifetime")
  O2_DEFINE_CONFIGURABLE(isDeepAngle, bool, true, "flag for applying deep angle")
  O2_DEFINE_CONFIGURABLE(cfgDeepAngle, float, 0.04, "deep angle cut")
  O2_DEFINE_CONFIGURABLE(removefaketrack, bool, true, "flag to remove fake kaon")
  O2_DEFINE_CONFIGURABLE(ConfFakeKaonCut, float, 0.1, "Cut based on track from momentum difference")
  O2_DEFINE_CONFIGURABLE(nsigmaCutTPC, float, 3, "nsigma tpc")
  O2_DEFINE_CONFIGURABLE(nsigmaCutTOF, float, 3, "nsigma tof")
  O2_DEFINE_CONFIGURABLE(cfgCutTOFBeta, float, 0.0, "TOF beta")
  O2_DEFINE_CONFIGURABLE(isTOFOnly, bool, false, "flag to select kaon with only TOF condition")

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

  using PIDTrack = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta>;
  using ResoV0s = aod::V0Datas;

  // std::unique_ptr<TFormula> sigmaFormula;

  void init(InitContext&)
  {
    // if (doprocessDataInvMass)
    // sigmaFormula = std::make_unique<TFormula>("sigmaFormula", cfgImSigmaFormula.value.c_str());
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
      if ((mcParticle.flagMcMatchGen() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) == 0 || mcParticle.daughtersIds().size() != 2) {
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
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isFakeKaon(T const& track)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (TMath::Abs(pglobal - ptpc) > ConfFakeKaonCut) {
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

    if (v0.mK0Short() < massK0Min || v0.mK0Short() > massK0Max) {
      return false;
    }
    if ((v0.qtarm() / std::abs(v0.alpha())) < qtArmenterosMinForK0) {
      return false;
    }
    if (v0.v0radius() > radiusMax || v0.v0radius() < radiusMin) {
      return false;
    }
    if (v0.v0cosPA() < cosPaMin) {
      return false;
    }
    if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
      return false;
    }
    if (v0.dcav0topv() > dcaV0ToPvMax) {
      return false;
    }
    if (std::abs(CtauK0s) > cMaxK0sLifeTime) {
      return false;
    }
    if (std::abs(v0.yK0Short()) > ConfV0Rap) {
      return false;
    }
    if (((std::abs(posTrack.tpcNSigmaPi()) > ConfDaughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > ConfDaughPIDCuts))) {
      return false;
    }
    if ((TMath::Abs(v0.dcapostopv()) < cMinV0DCAPi || TMath::Abs(v0.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0Cand>
  bool isSelectedV0AsLambda(Collision const& collision, const V0Cand& v0, int pid /*0: lambda, 1: antilambda*/)
  {
    const auto& posTrack = v0.template posTrack_as<PIDTrack>();
    const auto& negTrack = v0.template negTrack_as<PIDTrack>();

    float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda;

    if ((v0.mLambda() < massLambdaMin || v0.mLambda() > massLambdaMax) &&
        (v0.mAntiLambda() < massLambdaMin || v0.mAntiLambda() > massLambdaMax)) {
      return false;
    }
    if (v0.v0radius() > radiusMax || v0.v0radius() < radiusMin) {
      return false;
    }
    if (v0.v0cosPA() < cosPaMin) {
      return false;
    }
    if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
      return false;
    }
    if (v0.dcav0topv() > dcaV0ToPvMax) {
      return false;
    }
    if (pid == 0 && (TMath::Abs(v0.dcapostopv()) < cMinV0DCAPr || TMath::Abs(v0.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (TMath::Abs(v0.dcapostopv()) < cMinV0DCAPi || TMath::Abs(v0.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }
    if (pid == 0 && ((std::abs(posTrack.tpcNSigmaPr()) > ConfDaughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > ConfDaughPIDCuts))) {
      return false;
    }
    if (pid == 1 && ((std::abs(posTrack.tpcNSigmaPi()) > ConfDaughPIDCuts) || (std::abs(negTrack.tpcNSigmaPr()) > ConfDaughPIDCuts))) {
      return false;
    }
    if (std::abs(CtauLambda) > cMaxLambdaLifeTime) {
      return false;
    }
    if (std::abs(v0.yLambda()) > ConfV0Rap) {
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
    if (posTrack.tpcNClsCrossedRows() < tpcNClsCrossedRowsTrackMin || negTrack.tpcNClsCrossedRows() < tpcNClsCrossedRowsTrackMin) {
      return false;
    }
    if (posTrack.tpcCrossedRowsOverFindableCls() < 0.8 || negTrack.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }
    if (std::abs(v0.positiveeta()) > etaTrackMax || std::abs(v0.negativeeta()) > etaTrackMax) {
      return false;
    }
    if (v0.positivept() < ptTrackMin || v0.negativept() < ptTrackMin) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (cfgMomDepPID) {
      if (candidate.p() < 0.5) {
        if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
      } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
        return true;
      }
    } else if (!cfgMomDepPID) {
      if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      } else if (candidate.hasTOF() && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
        return true;
      }
    }
    return false;
  }

  // Generic 2-prong invariant mass method candidate finder. Only works for non-identical daughters of opposite charge for now.
  void processDataInvMass(aod::Collisions::iterator const& collision, aod::BCsWithTimestamps const&, aod::CFCollRefs const& cfcollisions, aod::CFTrackRefs const& cftracks, PIDTrack const& tracks, aod::V0Datas const& V0s)
  {
    if (cfcollisions.size() <= 0 || cftracks.size() <= 0)
      return; // rejected collision

    o2::aod::ITSResponse itsResponse;

    if (cfgDoPhi) {                           // Process Phi mesons
      for (const auto& cftrack1 : cftracks) { // Loop over first track
        const auto& p1 = tracks.iteratorAt(cftrack1.trackId() - tracks.begin().globalIndex());
        if (p1.sign() != 1) {
          continue;
        }
        if (ITSPIDSelection && p1.p() < ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) > -ITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p1) < ITSPIDNsigma.value)) { // Check ITS PID condition
          continue;
        }
        if (!selectionPID(p1)) {
          continue;
        }
        if (removefaketrack && isFakeKaon(p1)) { // Check if the track is a fake kaon
          continue;
        }
        for (const auto& cftrack2 : cftracks) {                 // Loop over second track
          if (cftrack2.globalIndex() == cftrack1.globalIndex()) // Skip if it's the same track as the first one
            continue;

          const auto& p2 = tracks.iteratorAt(cftrack2.trackId() - tracks.begin().globalIndex());
          if (p2.sign() != -1) {
            continue;
          }
          if (!selectionPID(p2)) {
            continue;
          }
          if (ITSPIDSelection && p2.p() < ITSPIDPthreshold.value && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) > -ITSPIDNsigma.value && itsResponse.nSigmaITS<o2::track::PID::Kaon>(p2) < ITSPIDNsigma.value)) { // Check ITS PID condition
            continue;
          }
          if (removefaketrack && isFakeKaon(p2)) { // Check if the track is a fake kaon
            continue;
          }
          if (!selectionPair(p1, p2)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector vec1(p1.pt(), p1.eta(), p1.phi(), cfgImPart1Mass);
          ROOT::Math::PtEtaPhiMVector vec2(p2.pt(), p2.eta(), p2.phi(), cfgImPart2Mass);
          ROOT::Math::PtEtaPhiMVector s = vec1 + vec2;
          if (s.pt() < cfgImCutPt || s.M() < cfgImMinInvMass || s.M() > cfgImMaxInvMass)
            continue;

          float phi = RecoDecay::constrainAngle(s.Phi(), 0.0f);
          output2ProngTracks(cfcollisions.begin().globalIndex(),
                             cftrack1.globalIndex(), cftrack2.globalIndex(), s.pt(), s.eta(), phi, s.M(), aod::cf2prongtrack::PhiToKK);
        } // end of loop over second track
      } // end of loop over first track
    } // end of processing Phi mesons

    if (cfgDoV0) {                    // Process V0 candidates (K0s, Lambdas, Anti-Lambdas)
      for (const auto& v0 : V0s) {    // Loop over V0 candidates
        if (!isV0TrackSelected(v0)) { // Quality selection for V0 prongs
          continue;
        }

        const auto& posTrack = v0.template posTrack_as<PIDTrack>();
        const auto& negTrack = v0.template negTrack_as<PIDTrack>();

        auto v0Type = 0;
        double massV0 = 0.0;
        if (isSelectedV0AsK0s(collision, v0)) { // candidate is K0s
          SETBIT(v0Type, aod::cf2prongtrack::K0stoPiPi);

          output2ProngTracks(cfcollisions.begin().globalIndex(),
                             posTrack.globalIndex(), negTrack.globalIndex(),
                             v0.pt(), v0.eta(), v0.phi(), v0.mK0Short(), v0Type);

          continue; // candidate is K0s, skip to next V0 candidate
        }

        bool LambdaTag = isSelectedV0AsLambda(collision, v0, 0);
        bool aLambdaTag = isSelectedV0AsLambda(collision, v0, 1);
        if (!LambdaTag && !aLambdaTag) { // neither Lambda nor Anti-Lambda
          continue;
        }
        // Note: candidate compatible with Lambda and Anti-Lambda hypothesis are counted twice (once for each hypothesis)
        if (LambdaTag) { // candidate is Lambda
          SETBIT(v0Type, aod::cf2prongtrack::LambdatoPPi);
          massV0 = v0.mLambda();
          output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                             v0.pt(), v0.eta(), v0.phi(), massV0, v0Type);
        } else if (aLambdaTag) { // candidate is Anti-lambda
          SETBIT(v0Type, aod::cf2prongtrack::AntiLambdatoPiP);
          massV0 = v0.mAntiLambda();
          output2ProngTracks(cfcollisions.begin().globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                             v0.pt(), v0.eta(), v0.phi(), massV0, v0Type);
        } // end of Lambda and Anti-Lambda processing
      } // end of loop over V0 candidates
    } // end of processing V0 candidates
  }
  PROCESS_SWITCH(Filter2Prong, processDataInvMass, "Process data generic 2-prong candidates with invariant mass method", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Filter2Prong>(cfgc, TaskName{"filter-2prong"})};
}
