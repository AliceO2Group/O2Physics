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

/// \file quarkGluonJetsProducer.cxx
/// \brief Produce a self-contained quark/gluon jet ML skim from JE derived data.

/// For every accepted detector-level charged jet the task writes:
///   * one jet row,
///   * one row per selected reconstructed jet constituent,
///   * low-level reconstructed PID information,
///   * MC-truth constituent information,
///   * leading and subleading quark/gluon candidates used for jet labelling.

/// \author Aleksandra Mulewicz <aleksandra.ewa.mulewicz@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace o2::aod
{
namespace qgmljet
{
DECLARE_SOA_COLUMN(FlavorLabel, flavorLabel, int32_t);
DECLARE_SOA_COLUMN(LeadingPartonPdg, leadingPartonPdg, int32_t);
DECLARE_SOA_COLUMN(LeadingPartonPt, leadingPartonPt, float);
DECLARE_SOA_COLUMN(LeadingPartonDeltaR, leadingPartonDeltaR, float);
DECLARE_SOA_COLUMN(SubleadingPartonPdg, subleadingPartonPdg, int32_t);
DECLARE_SOA_COLUMN(SubleadingPartonPt, subleadingPartonPt, float);
DECLARE_SOA_COLUMN(SubleadingPartonDeltaR, subleadingPartonDeltaR, float);
DECLARE_SOA_COLUMN(NPartonsInCone, nPartonsInCone, int32_t);
DECLARE_SOA_COLUMN(LabelAmbiguous, labelAmbiguous, bool);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetArea, jetArea, float);
DECLARE_SOA_COLUMN(JetRadius, jetRadius, float);
DECLARE_SOA_COLUMN(ZVertex, zVertex, float);
DECLARE_SOA_COLUMN(NConstituents, nConstituents, int32_t);
} // namespace qgmljet

DECLARE_SOA_TABLE(QGMLJets, "AOD", "QGMLJETS",
                  qgmljet::FlavorLabel,
                  qgmljet::LeadingPartonPdg,
                  qgmljet::LeadingPartonPt,
                  qgmljet::LeadingPartonDeltaR,
                  qgmljet::SubleadingPartonPdg,
                  qgmljet::SubleadingPartonPt,
                  qgmljet::SubleadingPartonDeltaR,
                  qgmljet::NPartonsInCone,
                  qgmljet::LabelAmbiguous,
                  qgmljet::JetPt,
                  qgmljet::JetEta,
                  qgmljet::JetPhi,
                  qgmljet::JetArea,
                  qgmljet::JetRadius,
                  qgmljet::ZVertex,
                  qgmljet::NConstituents);

using QGMLJet = QGMLJets::iterator;

namespace qgmlconst
{
DECLARE_SOA_INDEX_COLUMN(QGMLJet, qGMLJet);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaR, deltaR, float);
DECLARE_SOA_COLUMN(PtFraction, ptFraction, float);
DECLARE_SOA_COLUMN(Charge, charge, int32_t);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float);
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float);
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float);
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float);
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float);
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float);
DECLARE_SOA_COLUMN(TofBeta, tofBeta, float);
DECLARE_SOA_COLUMN(TruthPdg, truthPdg, int32_t);
DECLARE_SOA_COLUMN(TruthPt, truthPt, float);
DECLARE_SOA_COLUMN(TruthEta, truthEta, float);
DECLARE_SOA_COLUMN(TruthPhi, truthPhi, float);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
} // namespace qgmlconst

DECLARE_SOA_TABLE(QGMLConstituents, "AOD", "QGMLCONSTS",
                  qgmlconst::QGMLJetId,
                  qgmlconst::Pt,
                  qgmlconst::P,
                  qgmlconst::Eta,
                  qgmlconst::Phi,
                  qgmlconst::Px,
                  qgmlconst::Py,
                  qgmlconst::Pz,
                  qgmlconst::DeltaEta,
                  qgmlconst::DeltaPhi,
                  qgmlconst::DeltaR,
                  qgmlconst::PtFraction,
                  qgmlconst::Charge,
                  qgmlconst::HasTOF,
                  qgmlconst::DcaXY,
                  qgmlconst::DcaZ,
                  qgmlconst::TpcNSigmaPi,
                  qgmlconst::TpcNSigmaKa,
                  qgmlconst::TpcNSigmaPr,
                  qgmlconst::TofNSigmaPi,
                  qgmlconst::TofNSigmaKa,
                  qgmlconst::TofNSigmaPr,
                  qgmlconst::TofBeta,
                  qgmlconst::TruthPdg,
                  qgmlconst::TruthPt,
                  qgmlconst::TruthEta,
                  qgmlconst::TruthPhi,
                  qgmlconst::IsPhysicalPrimary);
} // namespace o2::aod

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;

using MCDJetEvents = aod::JetCollisionsMCD;
using ChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets,
                                 aod::ChargedMCDetectorLevelJetConstituents>;
using JetTrackRefs = soa::Join<aod::JetTracks, aod::JTrackPIs>;

using HadronTracksMC = soa::Join<aod::Tracks,
                                 aod::TracksExtra,
                                 aod::TrackSelection,
                                 aod::TrackSelectionExtension,
                                 aod::TracksDCA,
                                 aod::pidTPCFullPi,
                                 aod::pidTOFFullPi,
                                 aod::pidTPCFullKa,
                                 aod::pidTOFFullKa,
                                 aod::pidTPCFullPr,
                                 aod::pidTOFFullPr,
                                 aod::pidTOFbeta,
                                 aod::McTrackLabels>;

struct QuarkGluonJetsProducer {
  Produces<aod::QGMLJets> qgMLJets;
  Produces<aod::QGMLConstituents> qgMLConstituents;

  static constexpr double JetRadiusTolerance = 1.e-3;

  Configurable<std::string> eventSelections{"eventSelections", "selMCFull+NoITSROFrameBorder+IsGoodZvtxFT0vsPV", "JE derived-data event selections"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "skip minimum-bias gap subgenerator events"};
  Configurable<bool> applyRCTSelection{"applyRCTSelection", false, "apply RCT selection (normally false for MC training)"};
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT label used when applyRCTSelection=true"};
  std::vector<int> eventSelectionBits;

  Configurable<bool> isppRefAnalysis{"isppRefAnalysis", true, "pp reference jet acceptance"};
  Configurable<double> cfgEtaJetMax{"cfgEtaJetMax", 0.5, "maximum absolute jet eta in pp"};
  Configurable<double> minJetPt{"minJetPt", 10.0, "minimum detector-level jet pT"};
  Configurable<double> maxJetPt{"maxJetPt", 1000.0, "maximum detector-level jet pT"};
  Configurable<double> rJet{"rJet", 0.4, "selected jet radius"};
  Configurable<double> zVtx{"zVtx", 10.0, "maximum absolute collision z"};
  Configurable<bool> applyAreaCut{"applyAreaCut", false, "apply normalized jet-area cut"};
  Configurable<double> minNormalizedJetArea{"minNormalizedJetArea", 0.6, "minimum A/(pi R^2)"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from tracking edge"};

  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require PV-contributor track"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 70, "minimum TPC crossed rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi2/Ncl"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi2/Ncl"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi2/Ncl"};
  Configurable<double> minPt{"minPt", 0.15, "minimum selected-track pT"};
  Configurable<double> maxPt{"maxPt", 1000.0, "maximum selected-track pT"};
  Configurable<double> minEta{"minEta", -0.8, "minimum selected-track eta"};
  Configurable<double> maxEta{"maxEta", 0.8, "maximum selected-track eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.2, "maximum absolute DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.1, "maximum absolute DCAz"};

  Configurable<float> labelAmbiguityPtFraction{"labelAmbiguityPtFraction", 0.5f, "ambiguous if opposite-flavour subleading parton exceeds this leading-pT fraction"};
  Configurable<bool> requireGeneratorPartons{"requireGeneratorPartons", true, "use only MC partons produced by the event generator for Q/G labelling"};
  Configurable<int> partonGenStatus{"partonGenStatus", 0, "optional absolute generator-status filter for Q/G partons; 0 disables the status filter"};
  Configurable<bool> dropUnknownLabels{"dropUnknownLabels", false, "drop jets without quark/gluon label"};
  Configurable<bool> dropAmbiguousLabels{"dropAmbiguousLabels", false, "drop ambiguous labels"};

  Preslice<aod::JetParticles> mcParticlesPerCollision = aod::jmcparticle::mcCollisionId;

  enum JetFlavor { UnknownJet = 0,
                   QuarkJet = 1,
                   GluonJet = 2 };

  struct JetFlavorInfo {
    int flavor{UnknownJet};
    int leadingPartonPdg{0};
    float leadingPartonPt{-1.0f};
    float leadingPartonDeltaR{-1.0f};
    int subleadingPartonPdg{0};
    float subleadingPartonPt{-1.0f};
    float subleadingPartonDeltaR{-1.0f};
    int nPartonsInCone{0};
    bool ambiguous{false};
  };

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
  }

  template <typename TrackT>
  bool hasITSLayerHit(TrackT const& track, int layer) const
  {
    return TESTBIT(track.itsClusterMap(), layer - 1);
  }

  template <typename TrackT>
  bool passedTrackSelection(TrackT const& track)
  {
    if (requirePvContributor && !track.isPVContributor()) {
      return false;
    }
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (!hasITSLayerHit(track, 1) && !hasITSLayerHit(track, 2) && !hasITSLayerHit(track, 3)) {
      return false;
    }
    if (track.itsNCls() < minItsNclusters || track.tpcNClsCrossedRows() < minTpcNcrossedRows) {
      return false;
    }
    if (track.tpcChi2NCl() < minChiSquareTpc || track.tpcChi2NCl() > maxChiSquareTpc) {
      return false;
    }
    if (track.itsChi2NCl() > maxChiSquareIts) {
      return false;
    }
    if (track.eta() < minEta || track.eta() > maxEta || track.pt() < minPt || track.pt() > maxPt) {
      return false;
    }
    if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz) {
      return false;
    }
    return true;
  }

  template <typename JetT, typename ParticleRangeT>
  JetFlavorInfo getJetFlavorTag(JetT const& jet, ParticleRangeT const& particles)
  {
    JetFlavorInfo result;
    const double radius = static_cast<double>(jet.r()) / 100.0;

    for (auto const& particle : particles) {
      if (requireGeneratorPartons && !particle.producedByGenerator()) {
        continue;
      }
      if (partonGenStatus > 0 && std::abs(particle.getGenStatusCode()) != partonGenStatus) {
        continue;
      }

      const int absPdg = std::abs(particle.pdgCode());
      const bool isQuark = absPdg >= kDown && absPdg <= kTop;
      const bool isGluon = absPdg == kGluon;
      if (!isQuark && !isGluon) {
        continue;
      }

      const double dEta = particle.eta() - jet.eta();
      const double dPhi = RecoDecay::constrainAngle(particle.phi() - jet.phi(), -PI);
      const double dR = std::hypot(dEta, dPhi);
      if (dR >= radius) {
        continue;
      }

      ++result.nPartonsInCone;
      const float pt = particle.pt();
      if (pt > result.leadingPartonPt) {
        result.subleadingPartonPdg = result.leadingPartonPdg;
        result.subleadingPartonPt = result.leadingPartonPt;
        result.subleadingPartonDeltaR = result.leadingPartonDeltaR;
        result.leadingPartonPdg = particle.pdgCode();
        result.leadingPartonPt = pt;
        result.leadingPartonDeltaR = dR;
      } else if (pt > result.subleadingPartonPt) {
        result.subleadingPartonPdg = particle.pdgCode();
        result.subleadingPartonPt = pt;
        result.subleadingPartonDeltaR = dR;
      }
    }

    const int leadingAbs = std::abs(result.leadingPartonPdg);
    if (leadingAbs >= kDown && leadingAbs <= kTop) {
      result.flavor = QuarkJet;
    } else if (leadingAbs == kGluon) {
      result.flavor = GluonJet;
    }

    if (result.leadingPartonPt > 0.f && result.subleadingPartonPt > 0.f) {
      const int subAbs = std::abs(result.subleadingPartonPdg);
      const bool leadingQuark = leadingAbs >= kDown && leadingAbs <= kTop;
      const bool subQuark = subAbs >= kDown && subAbs <= kTop;
      result.ambiguous = (leadingQuark != subQuark) &&
                         result.subleadingPartonPt >= labelAmbiguityPtFraction * result.leadingPartonPt;
    }
    return result;
  }

  void process(MCDJetEvents::iterator const& collision,
               ChargedMCDJets const& detectorLevelJets,
               JetTrackRefs const&,
               HadronTracksMC const&,
               aod::JetParticles const& jetParticles,
               aod::McParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(
          collision,
          eventSelectionBits,
          static_cast<bool>(skipMBGapEvents),
          static_cast<bool>(applyRCTSelection),
          static_cast<std::string>(rctLabel))) {
      return;
    }
    if (std::abs(collision.posZ()) > zVtx) {
      return;
    }

    for (auto const& jet : detectorLevelJets) {
      const double jetRadius = static_cast<double>(jet.r()) / 100.0;
      if (std::abs(jetRadius - static_cast<double>(rJet)) > JetRadiusTolerance) {
        continue;
      }

      if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax) {
        continue;
      }
      if (!isppRefAnalysis && (std::abs(jet.eta()) + jetRadius > maxEta - deltaEtaEdge)) {
        continue;
      }
      if (jet.pt() < minJetPt || jet.pt() > maxJetPt) {
        continue;
      }
      const double normalizedArea = jet.area() / (PI * jetRadius * jetRadius);
      if (applyAreaCut && normalizedArea < minNormalizedJetArea) {
        continue;
      }

      const int jetMcCollisionId = jet.collision_as<MCDJetEvents>().mcCollisionId();

      JetFlavorInfo flavor;
      if (jetMcCollisionId >= 0) {
        auto collisionParticles = jetParticles.sliceBy(mcParticlesPerCollision, jetMcCollisionId);
        flavor = getJetFlavorTag(jet, collisionParticles);
      }
      if (dropUnknownLabels && flavor.flavor == UnknownJet) {
        continue;
      }
      if (dropAmbiguousLabels && flavor.ambiguous) {
        continue;
      }

      int32_t nConstituents = 0;

      for (auto const& jtrack : jet.tracks_as<JetTrackRefs>()) {
        auto track = jtrack.track_as<HadronTracksMC>();

        if (!passedTrackSelection(track)) {
          continue;
        }

        ++nConstituents;
      }
      if (nConstituents == 0) {
        continue;
      }

      qgMLJets(
        flavor.flavor,
        flavor.leadingPartonPdg,
        flavor.leadingPartonPt,
        flavor.leadingPartonDeltaR,
        flavor.subleadingPartonPdg,
        flavor.subleadingPartonPt,
        flavor.subleadingPartonDeltaR,
        flavor.nPartonsInCone,
        flavor.ambiguous,
        static_cast<float>(jet.pt()),
        static_cast<float>(jet.eta()),
        static_cast<float>(jet.phi()),
        static_cast<float>(jet.area()),
        static_cast<float>(jetRadius),
        static_cast<float>(collision.posZ()),
        nConstituents);

      const auto qgMLJetId = qgMLJets.lastIndex();

      for (auto const& jtrack : jet.tracks_as<JetTrackRefs>()) {
        auto track = jtrack.track_as<HadronTracksMC>();

        if (!passedTrackSelection(track)) {
          continue;
        }

        const double dEta = track.eta() - jet.eta();
        const double dPhi = RecoDecay::constrainAngle(track.phi() - jet.phi(), -PI);
        const double dR = std::hypot(dEta, dPhi);
        int32_t truthPdg = 0;
        float truthPt = std::numeric_limits<float>::quiet_NaN();
        float truthEta = std::numeric_limits<float>::quiet_NaN();
        float truthPhi = std::numeric_limits<float>::quiet_NaN();
        bool isPhysicalPrimary = false;

        if (track.has_mcParticle()) {
          auto truth = track.mcParticle_as<aod::McParticles>();
          truthPdg = truth.pdgCode();
          truthPt = truth.pt();
          truthEta = truth.eta();
          truthPhi = truth.phi();
          isPhysicalPrimary = truth.isPhysicalPrimary();
        }

        const float missingTOF = std::numeric_limits<float>::quiet_NaN();
        qgMLConstituents(
          qgMLJetId,
          static_cast<float>(track.pt()),
          static_cast<float>(track.p()),
          static_cast<float>(track.eta()),
          static_cast<float>(track.phi()),
          static_cast<float>(track.px()),
          static_cast<float>(track.py()),
          static_cast<float>(track.pz()),
          static_cast<float>(dEta),
          static_cast<float>(dPhi),
          static_cast<float>(dR),
          jet.pt() > 0.f ? static_cast<float>(track.pt() / jet.pt()) : 0.f,
          track.sign(),
          track.hasTOF(),
          static_cast<float>(track.dcaXY()),
          static_cast<float>(track.dcaZ()),
          static_cast<float>(track.tpcNSigmaPi()),
          static_cast<float>(track.tpcNSigmaKa()),
          static_cast<float>(track.tpcNSigmaPr()),
          track.hasTOF() ? static_cast<float>(track.tofNSigmaPi()) : missingTOF,
          track.hasTOF() ? static_cast<float>(track.tofNSigmaKa()) : missingTOF,
          track.hasTOF() ? static_cast<float>(track.tofNSigmaPr()) : missingTOF,
          track.hasTOF() ? static_cast<float>(track.beta()) : missingTOF,
          truthPdg,
          truthPt,
          truthEta,
          truthPhi,
          isPhysicalPrimary);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QuarkGluonJetsProducer>(cfgc)};
}
