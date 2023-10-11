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

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TString.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

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
} // namespace muon_source
DECLARE_SOA_TABLE(HfMuonSource, "AOD", "MUONSOURCE", muon_source::Pt, muon_source::DcaXY, muon_source::Source);
} // namespace o2::aod

struct HfTaskSingleMuonSource {
  Produces<aod::HfMuonSource> singleMuonSource;

  Configurable<bool> applyMcMask{"applyMcMask", true, "Flag of apply the mcMask selection"};
  Configurable<int> trackType{"trackType", 0, "Muon track type, validated values are 0, 1, 2, 3 and 4"};

  double etaLow = -3.6; // low edge of eta acceptance
  double etaUp = -2.5;  // up edge of eta acceptance
  double edgeZ = 10.0;  // edge of event position Z

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
      "SecondaryMu",
      "Hadron",
      "Unidentified"};

    AxisSpec axisDCA{5000, 0., 5., "DCA (cm)"};
    AxisSpec axisChi2{500, 0., 100., "#chi^{2} of MCH-MFT matching"};
    AxisSpec axisPt{200, 0., 100., "#it{p}_{T,reco} (GeV/#it{c})"};

    HistogramConfigSpec h3PtDCAChi2{HistType::kTH3F, {axisPt, axisDCA, axisChi2}};

    for (const auto& src : muonSources) {
      registry.add(Form("h3%sPtDCAChi2", src.Data()), "", h3PtDCAChi2);
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
      if (pdgAbs < 10)
        break; // Quark

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

      auto pdgData(TDatabasePDG::Instance()->GetParticle(mcPart.pdgCode()));
      if (pdgData && !pdgData->AntiParticle()) {
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
    return (isMuon(mask) && TESTBIT(mask, HasBeautyParent) && (!TESTBIT(mask, HasLightParent)));
  }

  // this muon comes directly from beauty decay
  bool isBeautyDecayMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && (!TESTBIT(mask, HasCharmParent)));
  }

  // this muon comes from non-prompt charm decay and does not have light flavor parent
  bool isNonpromptCharmMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && TESTBIT(mask, HasCharmParent));
  }

  // this muon comes from prompt charm decay and does not have light flavor parent
  bool isPromptCharmMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasLightParent)));
  }

  // this muon comes from light flavor quark decay
  bool isLightDecayMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasLightParent) && (!TESTBIT(mask, IsSecondary)));
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

    singleMuonSource(pt, dca, mask);

    if (isBeautyDecayMu(mask)) {
      registry.fill(HIST("h3BeautyDecayMuPtDCAChi2"), pt, dca, chi2);
    } else if (isNonpromptCharmMu(mask)) {
      registry.fill(HIST("h3NonpromptCharmMuPtDCAChi2"), pt, dca, chi2);
    } else if (isPromptCharmMu(mask)) {
      registry.fill(HIST("h3PromptCharmMuPtDCAChi2"), pt, dca, chi2);
    } else if (isLightDecayMu(mask)) {
      registry.fill(HIST("h3LightDecayMuPtDCAChi2"), pt, dca, chi2);
    } else if (isSecondaryMu(mask)) {
      registry.fill(HIST("h3SecondaryMuPtDCAChi2"), pt, dca, chi2);
    } else if (isHadron(mask)) {
      registry.fill(HIST("h3HadronPtDCAChi2"), pt, dca, chi2);
    } else if (isUnidentified(mask)) {
      registry.fill(HIST("h3UnidentifiedPtDCAChi2"), pt, dca, chi2);
    }
  }

  void process(MyCollisions::iterator const& collision,
               McMuons const& muons,
               aod::McParticles const&)
  {
    // event selections
    if (!collision.sel8()) {
      return;
    }
    if (std::abs(collision.posZ()) > edgeZ) {
      return;
    }

    for (const auto& muon : muons) {
      // muon selections
      if (muon.trackType() != trackType) {
        continue;
      }
      if (applyMcMask && (muon.mcMask() != 0)) {
        continue;
      }
      const auto eta(muon.eta());
      if ((eta >= etaUp) || (eta < etaLow)) {
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
