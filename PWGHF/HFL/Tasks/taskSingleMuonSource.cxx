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

#include "TMath.h"
#include "TString.h"
#include "TDatabasePDG.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
using McMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;

namespace
{
enum ParticleType {
  IsID = 0,            // this particle is identified
  IsMuon,              // this is a muon
  IsSecondary,         // this is a secondary particle
  HasLightParent,      // this particle has a light flavor parent
  HasCharmParent,      // this particle has a charm parent
  HasBeautyParent,     // this particle has a beauty parent
  HasQuarkoniumParent, // this particle has a quarkonium parent
  HasTauParent         // this particle has a tau parent
};

const auto nSrcs(7);
const TString muonSources[nSrcs]{
  "BeautyDecayMu",
  "NonpromptCharmMu",
  "PromptCharmMu",
  "LightDecayMu",
  "SecondaryMu",
  "Hadron",
  "Unidentified"};
} // namespace

struct HfTaskMuonSourceMc {
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
    AxisSpec axisDCA{500, 0., 5., "DCA (cm)"};
    AxisSpec axisChi2{500, 0., 100., "#chi^{2} of MCH-MFT matching"};
    AxisSpec axisPt{200, 0., 100., "#it{p}_{T,reco} (GeV/#it{c})"};

    HistogramConfigSpec h2PtDCA{HistType::kTH2F, {axisPt, axisDCA}};
    HistogramConfigSpec h2PtChi2{HistType::kTH2F, {axisPt, axisChi2}};

    for (auto const& src : muonSources) {
      registry.add(Form("h2%sPtDCA", src.Data()), "", h2PtDCA);
      registry.add(Form("h2%sPtChi2", src.Data()), "", h2PtChi2);
    }
  }

  // get the bitmask for muon sources identification
  uint8_t getMask(const McMuons::iterator& muon)
  {
    uint8_t mask(0);
    if (muon.has_mcParticle()) {
      SETBIT(mask, IsID);
    } else {
      return mask;
    }

    auto mcPart(muon.mcParticle());
    if (std::abs(mcPart.pdgCode()) == 13) {
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

      if (pdgAbs == 15) {
        // Tau
        SETBIT(mask, HasTauParent);
        continue;
      }

      const int pdgRem(pdgAbs % 100000);
      if ((pdgRem < 100) || (pdgRem >= 1000)) {
        continue;
      }
      // compute the flavor of constituent quark
      const int flv(pdgRem / TMath::Power(10, static_cast<int>(TMath::Log10(pdgRem))));
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
    return (TESTBIT(mask, IsID) && TESTBIT(mask, IsMuon));
  }

  // this muon is come from beauty decay and does not have light flavor parent
  bool isBeautyMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasBeautyParent) && (!TESTBIT(mask, HasLightParent)));
  }

  // this muon is directly come from beauty decay
  bool isBeautyDecayMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && (!TESTBIT(mask, HasCharmParent)));
  }

  // this muon is come from non-prompt charm decay and does not have light flavor parent
  bool isNonpromptCharmMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && TESTBIT(mask, HasCharmParent));
  }

  // this muon is come from prompt charm decay and does not have light flavor parent
  bool isPromptCharmMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasCharmParent) && (!TESTBIT(mask, HasBeautyParent)) && (!TESTBIT(mask, HasLightParent)));
  }

  // this muon is come from light flavor quark decay
  bool isLightDecayMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, HasLightParent) && (!TESTBIT(mask, IsSecondary)));
  }

  // this muon is come from transport
  bool isSecondaryMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, IsSecondary));
  }

  // this is a hadron
  bool isHadron(const uint8_t& mask)
  {
    return (TESTBIT(mask, IsID) && (!TESTBIT(mask, IsMuon)));
  }

  // this particle is undientified
  bool isUnidentified(const uint8_t& mask)
  {
    return (!TESTBIT(mask, IsID));
  }

  // fill the histograms of each particle types
  void fillHistograms(const McMuons::iterator& muon)
  {
    const auto mask(getMask(muon));
    const auto pt(muon.pt()), chi2(muon.chi2MatchMCHMFT());
    const auto dca(std::sqrt(std::pow(muon.fwdDcaX(), 2.) + std::pow(muon.fwdDcaY(), 2.)));

    if (isBeautyDecayMu(mask)) {
      registry.fill(HIST("h2BeautyDecayMuPtDCA"), pt, dca);
      registry.fill(HIST("h2BeautyDecayMuPtChi2"), pt, chi2);
    } else if (isNonpromptCharmMu(mask)) {
      registry.fill(HIST("h2NonpromptCharmMuPtDCA"), pt, dca);
      registry.fill(HIST("h2NonpromptCharmMuPtChi2"), pt, chi2);
    } else if (isPromptCharmMu(mask)) {
      registry.fill(HIST("h2PromptCharmMuPtDCA"), pt, dca);
      registry.fill(HIST("h2PromptCharmMuPtChi2"), pt, chi2);
    } else if (isLightDecayMu(mask)) {
      registry.fill(HIST("h2LightDecayMuPtDCA"), pt, dca);
      registry.fill(HIST("h2LightDecayMuPtChi2"), pt, chi2);
    } else if (isSecondaryMu(mask)) {
      registry.fill(HIST("h2SecondaryMuPtDCA"), pt, dca);
      registry.fill(HIST("h2SecondaryMuPtChi2"), pt, chi2);
    } else if (isHadron(mask)) {
      registry.fill(HIST("h2HadronPtDCA"), pt, dca);
      registry.fill(HIST("h2HadronPtChi2"), pt, chi2);
    } else if (isUnidentified(mask)) {
      registry.fill(HIST("h2UnidentifiedPtDCA"), pt, dca);
      registry.fill(HIST("h2UnidentifiedPtChi2"), pt, chi2);
    }
  }

  void process(MyCollisions::iterator const& collision, McMuons const& muons, aod::McParticles const&)
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
    adaptAnalysisTask<HfTaskMuonSourceMc>(cfgc),
  };
}
