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

/// \file dndetaBuilder.h
/// \brief configuration, selections and processing for charged-particle pseudorapidity densities on femto derived data
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
///
/// Needs a femto producer run in pass-through mode (collisions, tracks and, for mc, mc particles).
/// The event trigger (e.g. a femto pair or triplet) is passed in as a callable returning bool,
/// so the same builder is used by all dN/deta tasks independent of the trigger type.
///
/// Tracks are selected with the femto track bitmask via two partitions:
///   - global tracks (with TPC):   ITS + TPC quality bits (trackbuilder::ConfTrackSelectionDndetaWithTpc)
///   - ITS-only tracks (no TPC):   ITS quality bits only  (trackbuilder::ConfTrackSelectionDndetaWithoutTpc)
/// The QA of the selected tracks is done with one track histogram manager per track class.
/// This reproduces the global track selection, where the TPC cuts are only applied if the track has TPC.
/// A track is taken from the first partition only if it has TPC and from the second only if it has no TPC,
/// so no track is counted twice.

#ifndef PWGCF_FEMTO_CORE_DNDETABUILDER_H_
#define PWGCF_FEMTO_CORE_DNDETABUILDER_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/dndetaHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace o2::analysis::femto::dndetabuilder
{

/// general configuration
struct ConfDndeta : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Dndeta");
  o2::framework::Configurable<bool> requireInelGt0{"requireInelGt0", false, "Require INEL>0 (reco: at least one PV contributor in |eta|<1; gen: at least one charged primary in |eta|<1)"};
  o2::framework::Configurable<bool> doCorrelation{"doCorrelation", false, "Data only: fill multiplicity estimator correlations"};
  o2::framework::Configurable<bool> doMcEfficiency{"doMcEfficiency", true, "MC only: fill generated/reconstructed dN/deta for the efficiency"};
  o2::framework::Configurable<bool> doMcLoss{"doMcLoss", true, "MC only: fill event and signal loss histograms"};
  o2::framework::Configurable<bool> useGenTrigger{"useGenTrigger", false, "MC only: require the trigger also at generator level for the event/signal loss (only if the task supports it)"};
  o2::framework::Configurable<float> multEtaMax{"multEtaMax", 0.5f, "MC only: maximum |eta| of the generated multiplicity estimator used in the event/signal loss histograms"};
  // pT variation of the generated spectrum below the threshold, to estimate the uncertainty of the extrapolation to pT = 0
  // weight = p0 + p1 * pT for pT < ptVariationThreshold, 1 otherwise
  o2::framework::Configurable<float> ptVariationThreshold{"ptVariationThreshold", 0.1f, "MC only: generated particles below this pT (GeV/c) are reweighted for the pT variations"};
  o2::framework::Configurable<std::vector<double>> ptVariationUp{"ptVariationUp", {2., -10.}, "MC only: weight of the upward pT variation, p0 + p1 * pT"};
  o2::framework::Configurable<std::vector<double>> ptVariationDown{"ptVariationDown", {0.5, 5.}, "MC only: weight of the downward pT variation, p0 + p1 * pT"};
};

/// acceptance of the measurement, applied to reconstructed tracks and generated particles
/// the track quality is selected with the track partitions (bitmask)
struct ConfDndetaAcceptance : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("DndetaAcceptance");
  o2::framework::Configurable<float> etaMax{"etaMax", 1.0f, "Maximum |eta| (reconstructed tracks and generated particles). The eta range of the track partitions should be at least as large"};
  // regions with possible ITS acceptance problems, off by default
  o2::framework::Configurable<bool> applyExtraPhiCut{"applyExtraPhiCut", false, "Reject the phi regions below (reconstructed tracks and generated particles)"};
  o2::framework::Configurable<float> phiGapLow{"phiGapLow", 3.07666f, "Reject phiGapLow < phi < phiGapHigh"};
  o2::framework::Configurable<float> phiGapHigh{"phiGapHigh", 3.12661f, "Reject phiGapLow < phi < phiGapHigh"};
  o2::framework::Configurable<float> phiEdgeLow{"phiEdgeLow", 0.03f, "Reject phi <= phiEdgeLow"};
  o2::framework::Configurable<float> phiEdgeHigh{"phiEdgeHigh", 6.253f, "Reject phi >= phiEdgeHigh"};
};

constexpr double MinAbsChargeInUnitsOfThird = 3.; // TParticlePDG::Charge() is given in units of |e|/3
constexpr float InelGt0EtaMax = 1.f;              // INEL>0 is defined for |eta| < 1 (reco side: multNTracksPVeta1)
constexpr std::size_t NPtVariationParameters = 2;

class DndetaBuilder
{
 public:
  DndetaBuilder() = default;
  ~DndetaBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confDndeta,
            T2 const& confAcceptance,
            T3 const& confCollisionSelection,
            T4 const& confTracksGlobal,
            T5 const& confTrackBinningGlobal,
            T6 const& confTrackQaBinningGlobal,
            T7 const& confTracksItsOnly,
            T8 const& confTrackBinningItsOnly,
            T9 const& confTrackQaBinningItsOnly,
            T10 const& confBinning,
            bool isMc,
            bool doStrangeness)
  {
    // general
    mRequireInelGt0 = confDndeta.requireInelGt0.value;
    mDoCorrelation = confDndeta.doCorrelation.value;
    mDoMcEfficiency = confDndeta.doMcEfficiency.value;
    mDoMcLoss = confDndeta.doMcLoss.value;
    mUseGenTrigger = confDndeta.useGenTrigger.value;
    mMultEtaMax = confDndeta.multEtaMax.value;
    mPtVariationThreshold = confDndeta.ptVariationThreshold.value;
    mPtVariationUp = confDndeta.ptVariationUp.value;
    mPtVariationDown = confDndeta.ptVariationDown.value;
    if (mPtVariationUp.size() != NPtVariationParameters || mPtVariationDown.size() != NPtVariationParameters) {
      LOG(fatal) << "Dndeta.ptVariationUp and Dndeta.ptVariationDown need exactly " << NPtVariationParameters << " parameters (p0 + p1 * pT). Breaking...";
    }

    // collision selection
    mVtxZMin = confCollisionSelection.vtxZMin.value;
    mVtxZMax = confCollisionSelection.vtxZMax.value;
    mMultMin = confCollisionSelection.multMin.value;
    mMultMax = confCollisionSelection.multMax.value;
    mCentMin = confCollisionSelection.centMin.value;
    mCentMax = confCollisionSelection.centMax.value;
    mMagFieldMin = confCollisionSelection.magFieldMin.value;
    mMagFieldMax = confCollisionSelection.magFieldMax.value;
    mCollisionMask = confCollisionSelection.collisionMask.value;

    // acceptance
    mEtaMax = confAcceptance.etaMax.value;
    mApplyExtraPhiCut = confAcceptance.applyExtraPhiCut.value;
    mPhiGapLow = confAcceptance.phiGapLow.value;
    mPhiGapHigh = confAcceptance.phiGapHigh.value;
    mPhiEdgeLow = confAcceptance.phiEdgeLow.value;
    mPhiEdgeHigh = confAcceptance.phiEdgeHigh.value;

    // the dndeta track selections reuse the standard track selection object, check that it is configured for charged particle counting
    checkTrackSelection(confTracksGlobal);
    checkTrackSelection(confTracksItsOnly);

    // track QA (analysis + QA histograms, also in mc; the mc track classes are part of the dndeta histograms,
    // the mc histograms of the track histogram manager assume a single particle species)
    mTrackHistManagerGlobal.template init<modes::Mode::kReco_Qa>(registry, trackhistmanager::makeTrackQaHistSpecMap(confTrackBinningGlobal, confTrackQaBinningGlobal), confTracksGlobal, confTrackQaBinningGlobal);
    mTrackHistManagerItsOnly.template init<modes::Mode::kReco_Qa>(registry, trackhistmanager::makeTrackQaHistSpecMap(confTrackBinningItsOnly, confTrackQaBinningItsOnly), confTracksItsOnly, confTrackQaBinningItsOnly);

    // histograms
    uint32_t blocks = dndetahistmanager::kBlockEvent;
    if (isMc) {
      if (mDoMcEfficiency) {
        blocks |= dndetahistmanager::kBlockMcEfficiency;
      }
      if (mDoMcLoss) {
        blocks |= dndetahistmanager::kBlockLoss;
      }
    } else {
      blocks |= dndetahistmanager::kBlockData;
      if (mDoCorrelation) {
        blocks |= dndetahistmanager::kBlockCorrelation;
      }
      if (doStrangeness) {
        blocks |= dndetahistmanager::kBlockStrangeness;
      }
    }
    mHistManager.init(registry, dndetahistmanager::makeDndetaHistSpecMap(confBinning), blocks);
  }

  [[nodiscard]] bool useGenTrigger() const { return mUseGenTrigger; }

  // ---------------------------------------------------------------------------
  // data
  // ---------------------------------------------------------------------------

  /// \return true if the collision was selected (including the trigger)
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool processData(T1 const& col, T2& partitionWithTpc, T3& partitionWithoutTpc, T4& cache, T5 const& trigger)
  {
    if (!selectCollision(col, trigger)) {
      return false;
    }
    const float posZ = col.posZ();
    const float cent = col.cent();
    mHistManager.fillEvent(posZ, cent);

    auto tracksWithTpc = partitionWithTpc->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto tracksWithoutTpc = partitionWithoutTpc->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    int nch = 0;
    for (auto const& track : tracksWithTpc) {
      if (!isTrackAccepted(track, true)) {
        continue;
      }
      nch++;
      mTrackHistManagerGlobal.template fill<modes::Mode::kReco_Qa>(track, tracksWithTpc);
      mHistManager.fillDataTrack(posZ, cent, track.eta(), track.phi(), true);
    }
    for (auto const& track : tracksWithoutTpc) {
      if (!isTrackAccepted(track, false)) {
        continue;
      }
      nch++;
      mTrackHistManagerItsOnly.template fill<modes::Mode::kReco_Qa>(track, tracksWithoutTpc);
      mHistManager.fillDataTrack(posZ, cent, track.eta(), track.phi(), false);
    }
    if (mDoCorrelation) {
      mHistManager.fillCorrelation(static_cast<float>(nch), col.mult(), col.multFT0A(), col.multFT0C());
    }
    return true;
  }

  /// fill strangeness yields of a selected collision; lambdas and k0shorts are the already selected slices of this collision
  template <typename T1, typename T2, typename T3>
  void processStrangeness(T1 const& col, T2 const& lambdas, T3 const& k0shorts)
  {
    const float cent = col.cent();
    mHistManager.fillStrangenessCollision(col.posZ(), cent);
    for (auto const& lambda : lambdas) {
      mHistManager.fillLambda(cent, lambda.eta(), lambda.mass(), lambda.sign() < 0);
    }
    for (auto const& k0short : k0shorts) {
      mHistManager.fillK0short(cent, k0short.eta(), k0short.mass());
    }
  }

  // ---------------------------------------------------------------------------
  // mc
  // ---------------------------------------------------------------------------

  /// timeframe-level mc processing
  /// - reco: for every mc collision only the reconstructed collision with the most PV contributors is used
  ///         (emulates bestCollisionIndex); it has to pass the full event selection including the reco trigger
  /// - gen:  efficiency histograms for all generated collisions, loss histograms for the generated event class
  ///         (vertex z, optional INEL>0, optional generator-level trigger)
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
  void processMc(T1 const& cols,
                 T2 const& mcCols,
                 T3& partitionWithTpc,
                 T4& partitionWithoutTpc,
                 T5& cache,
                 T6 const& mcParticles,
                 T7 const& mcMothers,
                 T8 const& perMcColParticles,
                 T9& pdg,
                 T10 const& recoTrigger,
                 T11 const& genTrigger)
  {
    const auto nMcCols = static_cast<std::size_t>(mcCols.size());

    // best reconstructed collision per mc collision (largest number of PV contributors, as in mcCollisionExtra)
    std::vector<int64_t> bestReco(nMcCols, -1);
    std::vector<int> bestNumContrib(nMcCols, -1);
    for (auto const& col : cols) {
      if (!col.has_fMcCol()) {
        continue;
      }
      const auto mcColIndex = static_cast<std::size_t>(col.fMcColId());
      if (col.numContrib() > bestNumContrib[mcColIndex]) {
        bestNumContrib[mcColIndex] = col.numContrib();
        bestReco[mcColIndex] = col.globalIndex();
      }
    }

    // reconstructed side
    std::vector<int64_t> selectedReco(nMcCols, -1);
    for (auto const& col : cols) {
      if (!col.has_fMcCol() || bestReco[static_cast<std::size_t>(col.fMcColId())] != col.globalIndex()) {
        continue;
      }
      if (!selectCollision(col, [&]() {
            return recoTrigger(col);
          })) {
        continue;
      }
      selectedReco[static_cast<std::size_t>(col.fMcColId())] = col.globalIndex();
      if (mDoMcEfficiency) {
        auto tracksWithTpc = partitionWithTpc->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
        auto tracksWithoutTpc = partitionWithoutTpc->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
        fillRecoMc(col, tracksWithTpc, tracksWithoutTpc, mcParticles, mcMothers);
      }
    }

    // generated side
    for (auto const& mcCol : mcCols) {
      const auto mcColIndex = static_cast<std::size_t>(mcCol.globalIndex());
      const bool hasSelectedReco = selectedReco[mcColIndex] >= 0;
      auto particles = mcParticles.sliceBy(perMcColParticles, mcCol.globalIndex());

      // efficiency: centrality of the selected reconstructed collision (underflow if there is none)
      if (mDoMcEfficiency) {
        const float recoCent = hasSelectedReco ? cols.rawIteratorAt(selectedReco[mcColIndex]).cent() : -1.f;
        mHistManager.fillGenCollision(mcCol.posZ(), recoCent, hasSelectedReco);
        for (auto const& particle : particles) {
          if (!isGenParticleSelected(particle, pdg)) {
            continue;
          }
          mHistManager.fillGenParticle(mcCol.posZ(), recoCent, particle.eta(), particle.phi(), genType(particle.pdgCode()), hasSelectedReco,
                                       ptVariationWeight(particle.pt(), mPtVariationUp), ptVariationWeight(particle.pt(), mPtVariationDown));
        }
      }

      // event and signal loss: generated event class, true centrality
      if (!mDoMcLoss) {
        continue;
      }
      if (mcCol.posZ() < mVtxZMin || mcCol.posZ() > mVtxZMax) {
        continue;
      }
      int nChMult = 0;
      int nChEta1 = 0;
      for (auto const& particle : particles) {
        if (!isPhysicalPrimary(particle) || std::abs(particle.eta()) >= InelGt0EtaMax || !isCharged(particle.pdgCode(), pdg)) {
          continue;
        }
        nChEta1++;
        if (std::abs(particle.eta()) < mMultEtaMax) {
          nChMult++;
        }
      }
      if (mRequireInelGt0 && nChEta1 == 0) {
        continue;
      }
      if (mUseGenTrigger && !genTrigger(mcCol)) {
        continue;
      }
      const auto mult = static_cast<float>(nChMult);
      mHistManager.fillLossEvent(mcCol.cent(), mult, hasSelectedReco);
      for (auto const& particle : particles) {
        if (!isGenParticleSelected(particle, pdg)) {
          continue;
        }
        mHistManager.fillLossParticle(particle.eta(), mcCol.cent(), mult, hasSelectedReco);
      }
    }
  }

 private:
  /// warn about settings of the standard track selection object that do not make sense for charged particle counting
  template <typename T>
  static void checkTrackSelection(T const& conf)
  {
    const std::string name = conf.prefix;
    if (conf.chargeSign.value != 0) {
      LOG(warn) << name << ".chargeSign is " << conf.chargeSign.value << ", only one charge sign is counted. Set it to 0 for dN/deta.";
    }
    if (conf.maskLowMomentum.value != conf.maskHighMomentum.value) {
      LOG(warn) << name << ": maskLowMomentum and maskHighMomentum differ, the track selection depends on pidThres. Set both to the same track quality mask for dN/deta.";
    }
    if (conf.rejectionMaskLowMomentum.value != 0 || conf.rejectionMaskHighMomentum.value != 0) {
      LOG(warn) << name << ": rejection masks are set, tracks are rejected based on PID. Set them to 0 for dN/deta.";
    }
    LOG(info) << name << ": pT range " << conf.ptMin.value << " - " << conf.ptMax.value << " GeV/c, eta range " << conf.etaMin.value << " - " << conf.etaMax.value;
  }

  /// event selection with cutflow; the trigger is evaluated last, so the trigger QA
  /// only contains events that passed all other event selections
  template <typename T1, typename T2>
  bool selectCollision(T1 const& col, T2 const& trigger)
  {
    mHistManager.fillEventCutflow(dndetahistmanager::kEventAll);
    if (col.posZ() < mVtxZMin || col.posZ() > mVtxZMax) {
      return false;
    }
    mHistManager.fillEventCutflow(dndetahistmanager::kEventVtxZ);
    if (col.mult() < mMultMin || col.mult() > mMultMax ||
        col.cent() < mCentMin || col.cent() > mCentMax ||
        col.magField() < mMagFieldMin || col.magField() > mMagFieldMax ||
        (col.mask() & mCollisionMask) != mCollisionMask) {
      return false;
    }
    mHistManager.fillEventCutflow(dndetahistmanager::kEventCollisionSelection);
    if (mRequireInelGt0 && col.multNTracksPVeta1() <= 0) {
      return false;
    }
    mHistManager.fillEventCutflow(dndetahistmanager::kEventInelGt0);
    if (!trigger()) {
      return false;
    }
    mHistManager.fillEventCutflow(dndetahistmanager::kEventTrigger);
    return true;
  }

  /// weight of a pT variation, p[0] + p[1] * pT below the threshold, 1 above
  [[nodiscard]] double ptVariationWeight(float pt, std::vector<double> const& p) const
  {
    if (pt >= mPtVariationThreshold) {
      return 1.;
    }
    return p[0] + p[1] * pt;
  }

  [[nodiscard]] bool isInPhiGap(float phi) const
  {
    if (!mApplyExtraPhiCut) {
      return false;
    }
    return (phi > mPhiGapLow && phi < mPhiGapHigh) || phi <= mPhiEdgeLow || phi >= mPhiEdgeHigh;
  }

  /// track quality is already selected by the partition
  /// \param fromTpcPartition true if the track comes from the partition with TPC quality bits
  template <typename T>
  bool isTrackAccepted(T const& track, bool fromTpcPartition) const
  {
    // rows only stored to resolve a daughter index (e.g. V0 daughter from another collision) are no tracks of this collision
    if (track.isDaughterOnly()) {
      return false;
    }
    // tracks with TPC are only taken from the TPC partition and tracks without TPC only from the other one,
    // so the TPC cuts apply to every track with TPC and no track is counted twice
    if (track.hasTpc() != fromTpcPartition) {
      return false;
    }
    return std::abs(track.eta()) < mEtaMax && !isInPhiGap(track.phi());
  }

  template <typename T>
  static bool isPhysicalPrimary(T const& particle)
  {
    return particle.origin() == static_cast<datatypes::McOriginType>(modes::McOrigin::kPhysicalPrimary);
  }

  template <typename T>
  static bool isCharged(int pdgCode, T& pdg)
  {
    const auto* particle = pdg->GetParticle(pdgCode);
    return particle != nullptr && std::abs(particle->Charge()) >= MinAbsChargeInUnitsOfThird;
  }

  /// generated charged physical primary inside the acceptance
  template <typename T1, typename T2>
  bool isGenParticleSelected(T1 const& particle, T2& pdg) const
  {
    if (!isPhysicalPrimary(particle)) {
      return false;
    }
    if (std::abs(particle.eta()) >= mEtaMax || isInPhiGap(particle.phi())) {
      return false;
    }
    return isCharged(particle.pdgCode(), pdg);
  }

  static dndetahistmanager::GenType genType(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus:
        return dndetahistmanager::kGenPion;
      case PDG_t::kKPlus:
        return dndetahistmanager::kGenKaon;
      case PDG_t::kProton:
        return dndetahistmanager::kGenProton;
      default:
        return dndetahistmanager::kGenOther;
    }
  }

  static dndetahistmanager::RecoType recoTypeFromSpecies(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus:
        return dndetahistmanager::kRecoPion;
      case PDG_t::kKPlus:
        return dndetahistmanager::kRecoKaon;
      case PDG_t::kProton:
        return dndetahistmanager::kRecoProton;
      default:
        return dndetahistmanager::kRecoOther;
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillRecoMcTrack(T1 const& col, T2 const& track, T3 const& /*mcParticles*/, T4 const& /*mcMothers*/, std::vector<int64_t>& usedLabels)
  {
    const float posZ = col.posZ();
    const float cent = col.cent();
    if (!track.has_fMcParticle()) {
      mHistManager.fillRecoTrack(posZ, cent, track.eta(), track.phi(), dndetahistmanager::kRecoBkg);
      return;
    }
    auto mcParticle = track.template fMcParticle_as<T3>();
    // tracks matched to a particle from another mc collision are skipped (as in the legacy task)
    if (mcParticle.fMcColId() != col.fMcColId()) {
      return;
    }
    mHistManager.fillRecoTrack(posZ, cent, track.eta(), track.phi(), dndetahistmanager::kRecoAll);

    dndetahistmanager::RecoType type = dndetahistmanager::kRecoSecondary;
    if (isPhysicalPrimary(mcParticle)) {
      type = recoTypeFromSpecies(mcParticle.pdgCode());
    }
    if (mcParticle.has_fMcMother()) {
      const int motherPdg = mcParticle.template fMcMother_as<T4>().pdgCode();
      if (motherPdg == PDG_t::kK0Short || std::abs(motherPdg) == PDG_t::kLambda0) {
        type = dndetahistmanager::kRecoWeakDecay;
      }
    }
    if (std::find(usedLabels.begin(), usedLabels.end(), track.fMcParticleId()) != usedLabels.end()) {
      type = dndetahistmanager::kRecoFake;
    }
    usedLabels.push_back(track.fMcParticleId());
    mHistManager.fillRecoTrack(posZ, cent, mcParticle.eta(), mcParticle.phi(), type);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillRecoMc(T1 const& col, T2 const& tracksWithTpc, T3 const& tracksWithoutTpc, T4 const& mcParticles, T5 const& mcMothers)
  {
    mHistManager.fillRecoCollision(col.posZ(), col.cent());
    // labels are tracked over both partitions, so a fake is also found if the other track comes from the other partition
    std::vector<int64_t> usedLabels;
    for (auto const& track : tracksWithTpc) {
      if (isTrackAccepted(track, true)) {
        mTrackHistManagerGlobal.template fill<modes::Mode::kReco_Qa>(track, tracksWithTpc);
        fillRecoMcTrack(col, track, mcParticles, mcMothers, usedLabels);
      }
    }
    for (auto const& track : tracksWithoutTpc) {
      if (isTrackAccepted(track, false)) {
        mTrackHistManagerItsOnly.template fill<modes::Mode::kReco_Qa>(track, tracksWithoutTpc);
        fillRecoMcTrack(col, track, mcParticles, mcMothers, usedLabels);
      }
    }
  }

  dndetahistmanager::DndetaHistManager mHistManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixDndetaTrackGlobal> mTrackHistManagerGlobal;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixDndetaTrackItsOnly> mTrackHistManagerItsOnly;

  // general
  bool mRequireInelGt0 = false;
  bool mDoCorrelation = false;
  bool mDoMcEfficiency = true;
  bool mDoMcLoss = true;
  bool mUseGenTrigger = false;
  float mMultEtaMax = 0.f;
  float mPtVariationThreshold = 0.f;
  std::vector<double> mPtVariationUp;
  std::vector<double> mPtVariationDown;

  // collision selection
  float mVtxZMin = 0.f;
  float mVtxZMax = 0.f;
  float mMultMin = 0.f;
  float mMultMax = 0.f;
  float mCentMin = 0.f;
  float mCentMax = 0.f;
  int mMagFieldMin = 0;
  int mMagFieldMax = 0;
  datatypes::CollisionMaskType mCollisionMask = 0;

  // acceptance
  float mEtaMax = 0.f;
  bool mApplyExtraPhiCut = false;
  float mPhiGapLow = 0.f;
  float mPhiGapHigh = 0.f;
  float mPhiEdgeLow = 0.f;
  float mPhiEdgeHigh = 0.f;
};

} // namespace o2::analysis::femto::dndetabuilder

#endif // PWGCF_FEMTO_CORE_DNDETABUILDER_H_
