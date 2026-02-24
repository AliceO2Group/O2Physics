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
//
/// \file taskSingleMuonMultMc.cxx
/// \brief Task used to identify various sources of single muons and to calculated Axe in Monte Carlo simulation.
/// \author Md Samsul Islam <md.samsul.islam@cern.ch>, IITB

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include <Framework/ASoA.h>
#include <Framework/Configurable.h>

#include <TString.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::fwdtrack;

auto static constexpr MinCharge = 3.f;

namespace
{
enum ParticleType {
  IsIdentified = 0,    // this particle is identified
  IsMuon,              // this is a muon
  IsSecondary,         // this is a secondary particle
  HasTauParent,        // this particle has a tau parent
  HasWParent,          // this particle has a W parent
  HasZParent,          // this particle has a Z parent
  HasLightParent,      // this particle has a light flavor parent
  HasQuarkoniumParent, // this particle has a quarkonium parent
  HasCharmParent,      // this particle has a charm parent
  HasBeautyParent      // this particle has a beauty parent
};
} // namespace

struct HfTaskSingleMuonMultMc {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  o2::framework::Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<int> mcMaskSelection{"mcMaskSelection", 0, "McMask for correct match, valid values are 0 and 128"};
  Configurable<int> pdgQuark{"pdgQuark", 10, "pdg Max for Quarks"};
  Configurable<int> pdgRemMin{"pdgRemMin", 100, " Min. pdg Remnant for calculating Hadron pdg"};
  Configurable<int> pdgRemMax{"pdgRemMax", 10000, "Max. pdg Remnant for calculating Hadron pdg"};
  Configurable<int> flvMin{"flvMin", 4, "Min flavor of constituent quark"};
  Configurable<int> flvMax{"flvMax", 6, "Max flavor of constituent quark"};
  Configurable<int> flvLeading{"flvLeading", 10, "Base to extract leading flavor of constituent quark"};

  Configurable<int> charge{"charge", -1, "Muon track charge, validated values are 0, 1 and -1, 0 represents both 1 and -1"};
  Configurable<float> zVtxMax{"zVtxMax", 10., "maxium z of primary vertex [cm]"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.15, "minimum pt of tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "maximum pseudorapidity of tracks"};
  Configurable<float> etaMin{"etaMin", -3.6, "minimum pseudorapidity"};
  Configurable<float> etaMax{"etaMax", -2.5, "maximum pseudorapidity"};
  Configurable<float> pDcaMin{"pDcaMin", 324., "p*DCA value for small RAbsorb"};
  Configurable<float> pDcaMax{"pDcaMax", 594., "p*DCA value for large RAbsorb"};
  Configurable<float> rAbsorbMin{"rAbsorbMin", 17.6, "R at absorber end minimum value"};
  Configurable<float> rAbsorbMax{"rAbsorbMax", 89.5, "R at absorber end maximum value"};
  Configurable<float> rAbsorbMid{"rAbsorbMid", 26.5, "R at absorber end split point for different p*DCA selections"};

  using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using McMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
  using MyTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksIU, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;
  using McGenCollisions = soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos>;
  using McRecCollisions = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>>;

  // Filter Global Track for Multiplicty
  Filter trackFilter = ((nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin));

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(InitContext&)
  {
    const TString muonSources[]{
      "Identified",
      "Muon",
      "SecondaryMu",
      "LightDecayMu",
      "TauDecayMu",
      "WBosonDecayMu",
      "ZBosonDecayMu",
      "QuarkoniumDecayMu",
      "BeautyMu",
      "BeautyDecayMu",
      "NonpromptCharmMu",
      "PromptCharmMu",
      "OtherMu",
      "Hadron",
      "Unidentified"};

    AxisSpec axisColNumber{1, 0.5, 1.5, "Selected collisions"};
    AxisSpec axisMcMask{1001, -500.5, 500.5, "Mc Mask Selection"};
    AxisSpec axisMcLabel{1001, -500.5, 500.5, "McLabel"};
    AxisSpec axisNCh{500, 0.5, 500.5, "#it{N}_{ch}"};
    AxisSpec axisNtrk{500, 0.5, 500.5, "#it{N}_{trk}"};
    AxisSpec axisPt{5000, 0., 500., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{1000, -5., 5., "#it{#eta}"};
    AxisSpec axisRAbsorb{1000, 0., 100., "#it{R}_{Absorb} (cm)"};
    AxisSpec axisDCA{500, 0., 5., "#it{DCA}_{xy} (cm)"};
    AxisSpec axisChi2{1000, 0., 1000., "MCH-MFT matching #chi^{2}"};
    AxisSpec axisPDca{100000, 0, 100000, "#it{p} #times DCA (GeV/#it{c} * cm)"};
    AxisSpec axisDeltaPt{10000, -50, 50, "#Delta #it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisTrackType{8, -1.5, 6.5, "TrackType"};
    AxisSpec axisZ{1000, -50, 50, "V_{z} axis"};
    AxisSpec axisCharge{21, -10.5, 10.5, "charge"};
    AxisSpec axisPDG{201, -100.5, 100.5, "PDG"};
    AxisSpec axisMuonMask{15, -0.5, 14.5, "MuonMask"};

    HistogramConfigSpec hNEventGen{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec hNEventRec{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec hParticleGen{HistType::kTHnSparseF, {axisPt, axisEta}};
    HistogramConfigSpec hParticleRec{HistType::kTHnSparseF, {axisPt, axisEta}};
    HistogramConfigSpec hTrackResponse{HistType::kTH2F, {axisNCh, axisNtrk}};

    HistogramConfigSpec hMcMask{HistType::kTH1F, {axisMcMask}};
    HistogramConfigSpec hMuMcLabel{HistType::kTH1F, {axisMcLabel}};
    HistogramConfigSpec hMuTrackType{HistType::kTH1F, {axisTrackType}};
    HistogramConfigSpec hNEventGenMu{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec hNEventRecMu{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec hGenMuVar{HistType::kTHnSparseF, {axisPDG, axisPt, axisEta}};
    HistogramConfigSpec hRecMuVar{HistType::kTHnSparseF, {axisPDG, axisPt, axisEta}};
    HistogramConfigSpec hMuAfterAccCuts{HistType::kTHnSparseF, {axisPDG, axisPt, axisEta, axisRAbsorb, axisDCA, axisPDca, axisChi2, axisTrackType}};

    registry.add("hNEventRec", "", hNEventRec);
    registry.add("hNEventGen", "", hNEventGen);
    registry.add("hParticleGen", "Generated particles ", hParticleGen);
    registry.add("hParticleRec", "Reconstructed particles", hParticleRec);
    registry.add("hTrackResponse", "Generation-Recontsruction Response", hTrackResponse);

    registry.add("hMcMask", "Muon mcMask", hMcMask);
    registry.add("hMuMcLabel", "Muon mcLabel", hMuMcLabel);
    registry.add("hMuTrackType", "PDG code", hMuTrackType);
    registry.add("hNEventGenMu", "Muon Generated", hNEventGenMu);
    registry.add("hNEventRecMu", "Muon Reconstruced", hNEventRecMu);
    registry.add("hGenMuVar", "Muon Generated Variables", hGenMuVar);
    registry.add("hRecMuVar", "Muon Reconstructed Variables", hRecMuVar);
    registry.add("hMuAfterAccCuts", "Muon Reconstructed Variables", hMuAfterAccCuts);

    for (const auto& src : muonSources) {
      registry.add(Form("h2%s", src.Data()), "", {HistType::kTH2F, {axisPt, axisEta}});
    }
  }

  // get the bitmask for muon sources identification
  uint16_t getMask(const McMuons::iterator& muon)
  {
    uint16_t mask(0);
    if (!muon.has_mcParticle()) {
      return mask;
    }
    SETBIT(mask, IsIdentified);

    auto mcPart(muon.mcParticle());
    if (std::abs(mcPart.pdgCode()) != kMuonMinus) {
      return mask;
    }
    SETBIT(mask, IsMuon);

    while (mcPart.has_mothers()) {
      mcPart = *(mcPart.mothers_first_as<aod::McParticles>());
      const int pdgAbs = std::abs(mcPart.pdgCode());

      // Stop at quark
      if (pdgAbs < pdgQuark) {
        break;
      }

      // Produced in transport code
      if (!mcPart.producedByGenerator()) {
        SETBIT(mask, IsSecondary);
        continue;
      }

      // Tau parent
      if (pdgAbs == kTauMinus) {
        SETBIT(mask, HasTauParent);
        continue;
      }

      // W boson
      if (pdgAbs == kWPlus) {
        SETBIT(mask, HasWParent);
        continue;
      }

      // Z boson
      if (pdgAbs == kZ0) {
        SETBIT(mask, HasZParent);
        continue;
      }

      const int pdgRem(pdgAbs % 100000);

      if (pdgRem == kProton) {
        continue;
      } // Beam particle

      if ((pdgRem < pdgRemMin) || (pdgRem >= pdgRemMax)) {
        continue;
      }
      // compute the leding flavor of constituent quark
      int flv(pdgRem);
      while (flv >= flvLeading) {
        flv /= flvLeading;
      }

      if (flv > flvMax) {
        // no more than 6 flavors
        continue;
      }
      if (flv < flvMin) {
        // light flavor
        SETBIT(mask, HasLightParent);
        continue;
      }

      auto pdgData(pdg->GetParticle(mcPart.pdgCode()));
      if (pdgData && (pdgData->AntiParticle() == nullptr)) {
        SETBIT(mask, HasQuarkoniumParent);
        continue;
      } else if (flv == flvMin) {
        SETBIT(mask, HasCharmParent);
        continue;
      } else {
        SETBIT(mask, HasBeautyParent);
        continue;
      }
    }

    return mask;
  }

  // particle has an associated MC particle
  bool isIdentified(const uint16_t mask)
  {
    return (TESTBIT(mask, IsIdentified));
  }
  // this particle is muon
  bool isMuon(const uint16_t mask)
  {
    return (TESTBIT(mask, IsIdentified) && TESTBIT(mask, IsMuon));
  }

  // this muon comes from transport
  bool isSecondaryMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, IsSecondary));
  }

  // this muon comes from light flavor quark decay
  bool isLightDecayMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasLightParent) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from tau decays
  bool isTauDecayMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasTauParent) && (!TESTBIT(mask, HasWParent)) && (!TESTBIT(mask, HasZParent)) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasCharmParent)));
  }

  // this muon comes from W+- decay
  bool isWBosonDecayMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasWParent) && (!TESTBIT(mask, HasZParent)) && (!TESTBIT(mask, HasTauParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from Z decay
  bool isZBosonDecayMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasZParent) && (!TESTBIT(mask, HasWParent)) && (!TESTBIT(mask, HasTauParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from quarkonium decay
  bool isQuarkoniumDecayMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasQuarkoniumParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasCharmParent)) && (!TESTBIT(mask, HasLightParent)));
  }

  // this muon comes from beauty decay and does not have light flavor parent
  bool isBeautyMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasBeautyParent) && (!TESTBIT(mask, HasQuarkoniumParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes directly from beauty decay
  bool isBeautyDecayMu(const uint16_t mask)
  {
    return (isBeautyMu(mask) && (!TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasQuarkoniumParent))) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from non-prompt charm decay and does not have light flavor parent
  bool isNonpromptCharmMu(const uint16_t mask)
  {
    return (isBeautyMu(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasQuarkoniumParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from prompt charm decay and does not have light flavor parent
  bool isPromptCharmMu(const uint16_t mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasQuarkoniumParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon comes from other sources which have not classified above.
  bool isOtherMu(const uint16_t mask)
  {
    return (isMuon(mask) && (!isSecondaryMu(mask)) && (!isLightDecayMu(mask)) && (!isTauDecayMu(mask)) && (!isWBosonDecayMu(mask)) && (!isZBosonDecayMu(mask)) && (!isQuarkoniumDecayMu(mask)) && (!isBeautyMu(mask)) && (!isPromptCharmMu(mask)));
  }

  // this is a hadron
  bool isHadron(const uint16_t mask)
  {
    return (TESTBIT(mask, IsIdentified) && (!TESTBIT(mask, IsMuon)));
  }

  // this particle is unidentified
  bool isUnidentified(const uint16_t mask)
  {
    return ((!TESTBIT(mask, IsIdentified)));
  }

  // fill the histograms of each particle types
  void fillHistograms(const McMuons::iterator& muon)
  {
    const auto mask(getMask(muon));
    const auto pt(muon.pt()), eta(muon.eta());

    if (isIdentified(mask)) {
      registry.fill(HIST("h2Identified"), pt, eta);
    }
    if (isMuon(mask)) {
      registry.fill(HIST("h2Muon"), pt, eta);
    }
    if (isSecondaryMu(mask)) {
      registry.fill(HIST("h2SecondaryMu"), pt, eta);
    }
    if (isLightDecayMu(mask)) {
      registry.fill(HIST("h2LightDecayMu"), pt, eta);
    }
    if (isTauDecayMu(mask)) {
      registry.fill(HIST("h2TauDecayMu"), pt, eta);
    }
    if (isWBosonDecayMu(mask)) {
      registry.fill(HIST("h2WBosonDecayMu"), pt, eta);
    }
    if (isZBosonDecayMu(mask)) {
      registry.fill(HIST("h2ZBosonDecayMu"), pt, eta);
    }
    if (isQuarkoniumDecayMu(mask)) {
      registry.fill(HIST("h2QuarkoniumDecayMu"), pt, eta);
    }
    if (isBeautyMu(mask)) {
      registry.fill(HIST("h2BeautyMu"), pt, eta);
    }
    if (isBeautyDecayMu(mask)) {
      registry.fill(HIST("h2BeautyDecayMu"), pt, eta);
    }
    if (isNonpromptCharmMu(mask)) {
      registry.fill(HIST("h2NonpromptCharmMu"), pt, eta);
    }
    if (isPromptCharmMu(mask)) {
      registry.fill(HIST("h2PromptCharmMu"), pt, eta);
    }
    if (isOtherMu(mask)) {
      registry.fill(HIST("h2OtherMu"), pt, eta);
    }
    if (isHadron(mask)) {
      registry.fill(HIST("h2Hadron"), pt, eta);
    }
    if (isUnidentified(mask)) {
      registry.fill(HIST("h2Unidentified"), pt, eta);
    }
  }

  void process(McGenCollisions::iterator const& mccollision,
               McMuons const& muons,
               aod::McParticles const&,
               McRecCollisions const& collisions)
  {

    // event selections
    if (std::abs(mccollision.posZ()) > zVtxMax) {
      return;
    }

    registry.fill(HIST("hNEventGenMu"), 1);

    for (const auto& muon : muons) {
      if (!(muon.has_mcParticle())) {
        continue;
      }
      auto mcPart(muon.mcParticle());
      auto pdgGen(mcPart.pdgCode());
      auto etaGen(mcPart.eta());

      if (!(std::abs(pdgGen) == kMuonMinus)) {
        continue;
      }
      if ((etaGen >= etaMax) || (etaGen < etaMin)) {
        continue;
      }
      registry.fill(HIST("hGenMuVar"), pdgGen, mcPart.pt(), etaGen);
    }

    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > zVtxMax) {
        continue;
      }
      registry.fill(HIST("hNEventRecMu"), 1);

      for (const auto& muon : muons) {
        // muon selections
        registry.fill(HIST("hMcMask"), muon.mcMask());

        if (muon.mcMask() != mcMaskSelection) {
          continue;
        }

        if (!(muon.has_mcParticle())) {
          continue;
        }
        const auto pt(muon.pt()), eta(muon.eta()), rAbsorb(muon.rAtAbsorberEnd()), pDca(muon.pDca()), chi2(muon.chi2MatchMCHMFT());
        const auto dcaXY{RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY())};
        const auto muTrackType{muon.trackType()};
        const auto mcLabel(muon.mcParticleId());

        if (mcLabel < 0) {
          continue;
        }

        registry.fill(HIST("hMuMcLabel"), mcLabel);
        registry.fill(HIST("hMuTrackType"), muTrackType);

        auto mcparticle(muon.mcParticle());
        const auto pdg(mcparticle.pdgCode());

        if ((eta >= etaMax) || (eta < etaMin)) {
          continue;
        }
        if ((rAbsorb >= rAbsorbMax) || (rAbsorb < rAbsorbMin)) {
          continue;
        }
        if (pDca >= pDcaMax) {
          continue;
        }
        registry.fill(HIST("hRecMuVar"), pdg, pt, eta);
        registry.fill(HIST("hMuAfterAccCuts"), pdg, pt, eta, rAbsorb, dcaXY, pDca, chi2, muTrackType);

        fillHistograms(muon);
      }
    }
  }

  void processResTrack(McGenCollisions::iterator const& mccollision,
                       McRecCollisions const& collisions,
                       aod::McParticles const& particles,
                       MyTracks const& tracks)
  {
    // event selections
    if (std::abs(mccollision.posZ()) > zVtxMax) {
      return;
    }
    registry.fill(HIST("hNEventGen"), 1.);
    auto nP = 0;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (!particle.producedByGenerator()) {
        continue;
      }

      auto charge = 0.;
      auto* p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }

      if (std::abs(charge) < MinCharge) {
        continue;
      }
      if (particle.pt() < ptTrackMin || std::abs(particle.eta()) >= etaTrackMax) {
        continue;
      }

      registry.fill(HIST("hParticleGen"), particle.pt(), particle.eta());
      nP++;
    }
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > zVtxMax) {
        continue;
      }
      registry.fill(HIST("hNEventRec"), 1.);
      auto nTrk = 0;
      auto tracksample = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracksample) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        registry.fill(HIST("hParticleRec"), track.pt(), track.eta());
        ++nTrk;
      }
      if (nTrk < 1) {
        continue;
      }
      registry.fill(HIST("hTrackResponse"), nP, nTrk);
    }
  }
  PROCESS_SWITCH(HfTaskSingleMuonMultMc, processResTrack, "Process Track Reconstruction/Generation", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleMuonMultMc>(cfgc),
  };
}
