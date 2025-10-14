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
// \file taskSingleMuonSource.cxx
// \brief Task used to seperate single muons source in Monte Carlo simulation.
// \author Maolin Zhang <maolin.zhang@cern.ch>, CCNU

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TString.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
using McMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

namespace
{
enum ParticleType {
  IsIdentified = 0,    // this particle is identified
  IsMuon,              // this is a muon
  IsSecondary,         // this is a secondary particle
  HasLightParent,      // this particle has a light flavor parent
  HasCharmParent,      // this particle has a charm parent
  HasBeautyParent,     // this particle has a beauty parent
  HasQuarkoniumParent, // this particle has a quarkonium parent
  HasTauParent         // this particle has a tau parent
};
} // namespace

namespace o2::aod
{
namespace muon_source
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(Source, source, uint8_t);
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, float);
} // namespace muon_source
DECLARE_SOA_TABLE(HfMuonSource, "AOD", "MUONSOURCE", muon_source::Pt, muon_source::DcaXY, muon_source::Source, muon_source::DeltaPt);
} // namespace o2::aod

struct HfTaskSingleMuonSource {
  Produces<aod::HfMuonSource> singleMuonSource;

  Configurable<int> mcMaskSelection{"mcMaskSelection", 0, "McMask for correct match, valid values are 0 and 128"};
  Configurable<int> trackType{"trackType", 0, "Muon track type, validated values are 0, 1, 2, 3 and 4"};
  Configurable<int> charge{"charge", -1, "Muon track charge, validated values are 0, 1 and -1, 0 represents both 1 and -1"};

  double pDcaMax = 594.0; // p*DCA maximum value for large Rabs
  double rAbsMin = 26.5;  // R at absorber end minimum value
  double rAbsMax = 89.5;  // R at absorber end maximum value
  double etaLow = -3.6;   // low edge of eta acceptance
  double etaUp = -2.5;    // up edge of eta acceptance
  double edgeZ = 10.0;    // edge of event position Z

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(InitContext&)
  {
    const TString muonSources[]{
      "BeautyDecayMu",
      "NonpromptCharmMu",
      "PromptCharmMu",
      "LightDecayMu",
      "QuarkoniumDecayMu",
      "SecondaryMu",
      "Hadron",
      "Unidentified"};

    AxisSpec const axisColNumber{1, 0.5, 1.5, "Selected collisions"};
    AxisSpec const axisDCA{5000, 0., 5., "DCA (cm)"};
    AxisSpec const axisChi2{500, 0., 100., "#chi^{2} of MCH-MFT matching"};
    AxisSpec const axisPt{200, 0., 100., "#it{p}_{T,reco} (GeV/#it{c})"};
    AxisSpec const axisDeltaPt{1000, -50., 50., "#Delta #it{p}_{T} (GeV/#it{c})"};

    HistogramConfigSpec const h1ColNumber{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec const h1Pt{HistType::kTH1F, {axisPt}};
    HistogramConfigSpec const h2PtDCA{HistType::kTH2F, {axisPt, axisDCA}};
    HistogramConfigSpec const h2PtChi2{HistType::kTH2F, {axisPt, axisChi2}};
    HistogramConfigSpec const h2PtDeltaPt{HistType::kTH2F, {axisPt, axisDeltaPt}};

    registry.add("h1ColNumber", "", h1ColNumber);
    for (const auto& src : muonSources) {
      registry.add(Form("h1%sPt", src.Data()), "", h1Pt);
      registry.add(Form("h2%sPtDCA", src.Data()), "", h2PtDCA);
      registry.add(Form("h2%sPtChi2", src.Data()), "", h2PtChi2);
      registry.add(Form("h2%sPtDeltaPt", src.Data()), "", h2PtDeltaPt);
    }
  }

  // get the bitmask for muon sources identification
  uint8_t getMask(const McMuons::iterator& muon)
  {
    uint8_t mask(0);
    if (muon.has_mcParticle()) {
      SETBIT(mask, IsIdentified);
    } else {
      return mask;
    }

    auto mcPart(muon.mcParticle());
    if (std::abs(mcPart.pdgCode()) == kMuonMinus) {
      // Muon
      SETBIT(mask, IsMuon);
    } else {
      return mask;
    }

    while (mcPart.has_mothers()) {
      mcPart = *(mcPart.mothers_first_as<aod::McParticles>());

      const auto pdgAbs(std::abs(mcPart.pdgCode()));
      if (pdgAbs < 10) {
        break; // Quark
      }

      if (!mcPart.producedByGenerator()) { // Produced in transport code
        SETBIT(mask, IsSecondary);
        continue;
      }

      if (pdgAbs == kTauMinus) {
        // Tau
        SETBIT(mask, HasTauParent);
        continue;
      }

      const int pdgRem(pdgAbs % 100000);

      if (pdgRem == kProton) {
        continue;
      } // Beam particle

      if ((pdgRem < 100) || (pdgRem >= 10000)) {
        continue;
      }
      // compute the flavor of constituent quark
      const int flv(pdgRem / std::pow(10, static_cast<int>(std::log10(pdgRem))));
      if (flv > 6) {
        // no more than 6 flavors
        continue;
      }
      if (flv < 4) {
        // light flavor
        SETBIT(mask, HasLightParent);
        continue;
      }

      auto* pdgData(TDatabasePDG::Instance()->GetParticle(mcPart.pdgCode()));
      if ((pdgData != nullptr) && (pdgData->AntiParticle() == nullptr)) {
        SETBIT(mask, HasQuarkoniumParent);
      } else if (flv == 4) {
        SETBIT(mask, HasCharmParent);
      } else {
        SETBIT(mask, HasBeautyParent);
      }
    }

    return mask;
  }

  // this particle is muon
  bool isMuon(const uint8_t& mask)
  {
    return (TESTBIT(mask, IsIdentified) && TESTBIT(mask, IsMuon));
  }

  // this muon comes from beauty decay and does not have light flavor parent
  bool isBeautyMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasBeautyParent) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, HasQuarkoniumParent)));
  }

  // this muon comes directly from beauty decay
  bool isBeautyDecayMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && (!TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasQuarkoniumParent))));
  }

  // this muon comes from non-prompt charm decay and does not have light flavor parent
  bool isNonpromptCharmMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasQuarkoniumParent)));
  }

  // this muon comes from prompt charm decay and does not have light flavor parent
  bool isPromptCharmMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasLightParent)) && (!TESTBIT(mask, HasQuarkoniumParent)));
  }

  // this muon comes from light flavor quark decay
  bool isLightDecayMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasLightParent) && (!TESTBIT(mask, IsSecondary)) && (!TESTBIT(mask, HasQuarkoniumParent)));
  }

  // this muon comes from quarkonium decay
  bool isQuarkoniumDecayMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasQuarkoniumParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasCharmParent)));
  }

  // this muon comes from transport
  bool isSecondaryMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, IsSecondary));
  }

  // this is a hadron
  bool isHadron(const uint8_t& mask)
  {
    return (TESTBIT(mask, IsIdentified) && (!TESTBIT(mask, IsMuon)));
  }

  // this particle is unidentified
  bool isUnidentified(const uint8_t& mask)
  {
    return (!TESTBIT(mask, IsIdentified));
  }

  // fill the histograms of each particle types
  void fillHistograms(const McMuons::iterator& muon)
  {
    const auto mask(getMask(muon));
    const auto pt(muon.pt()), chi2(muon.chi2MatchMCHMFT());
    const auto dca(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));

    if (trackType == 0 || trackType == 2) {
      if (!muon.has_matchMCHTrack()) {
        return;
      }
      const auto muonType3 = muon.matchMCHTrack_as<McMuons>();
      const auto deltaPt = muonType3.pt() - pt;

      singleMuonSource(pt, dca, mask, deltaPt);

      if (isBeautyDecayMu(mask)) {
        registry.fill(HIST("h2BeautyDecayMuPtDCA"), pt, dca);
        registry.fill(HIST("h2BeautyDecayMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2BeautyDecayMuPtDeltaPt"), pt, deltaPt);
      } else if (isNonpromptCharmMu(mask)) {
        registry.fill(HIST("h2NonpromptCharmMuPtDCA"), pt, dca);
        registry.fill(HIST("h2NonpromptCharmMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2NonpromptCharmMuPtDeltaPt"), pt, deltaPt);
      } else if (isPromptCharmMu(mask)) {
        registry.fill(HIST("h2PromptCharmMuPtDCA"), pt, dca);
        registry.fill(HIST("h2PromptCharmMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2PromptCharmMuPtDeltaPt"), pt, deltaPt);
      } else if (isLightDecayMu(mask)) {
        registry.fill(HIST("h2LightDecayMuPtDCA"), pt, dca);
        registry.fill(HIST("h2LightDecayMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2LightDecayMuPtDeltaPt"), pt, deltaPt);
      } else if (isQuarkoniumDecayMu(mask)) {
        registry.fill(HIST("h2QuarkoniumDecayMuPtDCA"), pt, dca);
        registry.fill(HIST("h2QuarkoniumDecayMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2QuarkoniumDecayMuPtDeltaPt"), pt, deltaPt);
      } else if (isSecondaryMu(mask)) {
        registry.fill(HIST("h2SecondaryMuPtDCA"), pt, dca);
        registry.fill(HIST("h2SecondaryMuPtChi2"), pt, chi2);
        registry.fill(HIST("h2SecondaryMuPtDeltaPt"), pt, deltaPt);
      } else if (isHadron(mask)) {
        registry.fill(HIST("h2HadronPtDCA"), pt, dca);
        registry.fill(HIST("h2HadronPtChi2"), pt, chi2);
        registry.fill(HIST("h2HadronPtDeltaPt"), pt, deltaPt);
      } else if (isUnidentified(mask)) {
        registry.fill(HIST("h2UnidentifiedPtDCA"), pt, dca);
        registry.fill(HIST("h2UnidentifiedPtChi2"), pt, chi2);
        registry.fill(HIST("h2UnidentifiedPtDeltaPt"), pt, deltaPt);
      }
    } else {
      if (isBeautyDecayMu(mask)) {
        registry.fill(HIST("h1BeautyDecayMuPt"), pt);
      } else if (isNonpromptCharmMu(mask)) {
        registry.fill(HIST("h1NonpromptCharmMuPt"), pt);
      } else if (isPromptCharmMu(mask)) {
        registry.fill(HIST("h1PromptCharmMuPt"), pt);
      } else if (isLightDecayMu(mask)) {
        registry.fill(HIST("h1LightDecayMuPt"), pt);
      } else if (isQuarkoniumDecayMu(mask)) {
        registry.fill(HIST("h1QuarkoniumDecayMuPt"), pt);
      } else if (isSecondaryMu(mask)) {
        registry.fill(HIST("h1SecondaryMuPt"), pt);
      } else if (isHadron(mask)) {
        registry.fill(HIST("h1HadronPt"), pt);
      } else if (isUnidentified(mask)) {
        registry.fill(HIST("h1UnidentifiedPt"), pt);
      }
    }
  }

  void process(MyCollisions::iterator const& collision,
               McMuons const& muons,
               aod::McParticles const&)
  {
    // event selections
    if (std::abs(collision.posZ()) > edgeZ) {
      return;
    }
    registry.fill(HIST("h1ColNumber"), 1.);

    for (const auto& muon : muons) {
      // muon selections
      if (muon.trackType() != trackType) {
        continue;
      }
      if (trackType == 0 && muon.mcMask() != mcMaskSelection) {
        continue;
      }
      const auto eta(muon.eta()), pDca(muon.pDca()), rAbs(muon.rAtAbsorberEnd());
      if ((eta >= etaUp) || (eta < etaLow)) {
        continue;
      }
      if ((rAbs >= rAbsMax) || (rAbs < rAbsMin)) {
        continue;
      }
      if (pDca >= pDcaMax) {
        continue;
      }
      if ((muon.chi2() >= 1e6) || (muon.chi2() < 0)) {
        continue;
      }
      if (charge != 0 && muon.sign() != charge) {
        continue;
      }
      fillHistograms(muon);
    } // loop over muons
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleMuonSource>(cfgc),
  };
}
