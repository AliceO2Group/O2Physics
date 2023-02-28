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
#include "ReconstructionDataFormats/TrackFwd.h"
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
  ParticleIsID = 0,
  ParticleIsMuon,
  ParticleIsSecondary,
  ParticleHasLightParent,
  ParticleHasCharmParent,
  ParticleHasBeautyParent,
  ParticleHasQuarkoniumParent,
  ParticleHasTauParent
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

struct HfTaskMuonSourceMC {
  Configurable<bool> mApplyMcMask{"mApplyMcMask", true, "Flag of apply the mcMask selection"};
  Configurable<int> mTrackType{"mTrackType", 0, "Muon track type, validated values are 0, 1, 2, 3 and 4"};

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

    for (auto& src : muonSources) {
      registry.add(Form("h2%sPtDCA", src.Data()), "", h2PtDCA);
      registry.add(Form("h2%sPtChi2", src.Data()), "", h2PtChi2);
    }
  }

  uint8_t getMask(const McMuons::iterator& muon)
  {
    uint8_t mask(0);
    if (muon.has_mcParticle()) {
      SETBIT(mask, ParticleIsID);
    } else {
      return mask;
    }

    auto mcPart(muon.mcParticle());
    if (abs(mcPart.pdgCode()) == 13) {
      SETBIT(mask, ParticleIsMuon);
    } else {
      return mask;
    }

    while (mcPart.has_mothers()) {
      mcPart = *(mcPart.mothers_first_as<aod::McParticles>());

      const auto pdgAbs(abs(mcPart.pdgCode()));
      if (pdgAbs < 10)
        break; // Quark

      if (!mcPart.producedByGenerator()) { // Produced in transport code
        SETBIT(mask, ParticleIsSecondary);
        continue;
      }

      if (pdgAbs == 15) {
        SETBIT(mask, ParticleHasTauParent);
        continue;
      }

      const int pdgRem(pdgAbs % 100000);
      if ((pdgRem < 100) || (pdgRem >= 1000))
        continue;
      const int flv(pdgRem / TMath::Power(10, int(TMath::Log10(pdgRem))));
      if (flv > 6) {
        continue;
      }
      if (flv < 4) {
        SETBIT(mask, ParticleHasLightParent);
        continue;
      }

      auto pdgData(TDatabasePDG::Instance()->GetParticle(mcPart.pdgCode()));
      if (pdgData && !pdgData->AntiParticle()) {
        SETBIT(mask, ParticleHasQuarkoniumParent);
      } else if (flv == 4) {
        SETBIT(mask, ParticleHasCharmParent);
      } else {
        SETBIT(mask, ParticleHasBeautyParent);
      }
    }

    return mask;
  }

  bool isMuon(const uint8_t& mask)
  {
    return (TESTBIT(mask, ParticleIsID) && TESTBIT(mask, ParticleIsMuon));
  }

  bool isBeautyMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, ParticleHasBeautyParent) && (!TESTBIT(mask, ParticleHasLightParent)));
  }

  bool isBeautyDecayMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && (!TESTBIT(mask, ParticleHasCharmParent)));
  }

  bool isNonpromptCharmMu(const uint8_t& mask)
  {
    return (isBeautyMu(mask) && TESTBIT(mask, ParticleHasCharmParent));
  }

  bool isPromptCharmMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, ParticleHasCharmParent) && (!TESTBIT(mask, ParticleHasBeautyParent)) && (!TESTBIT(mask, ParticleHasLightParent)));
  }

  bool isLightDecayMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, ParticleHasLightParent) && (!TESTBIT(mask, ParticleIsSecondary)));
  }

  bool isSecondaryMu(const uint8_t& mask)
  {
    return (isMuon(mask) && TESTBIT(mask, ParticleIsSecondary));
  }

  bool isHadron(const uint8_t& mask)
  {
    return (TESTBIT(mask, ParticleIsID) && (!TESTBIT(mask, ParticleIsMuon)));
  }

  bool isUnidentified(const uint8_t& mask)
  {
    return (!TESTBIT(mask, ParticleIsID));
  }

  void fillHistograms(const McMuons::iterator& muon)
  {
    const auto mask(getMask(muon));
    const auto pt(muon.pt()), chi2(muon.chi2MatchMCHMFT());
    const auto dca(std::sqrt(pow(muon.fwdDcaX(), 2.) + pow(muon.fwdDcaY(), 2.)));

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
    if (!collision.sel8()) {
      return;
    }
    if (abs(collision.posZ()) > 10.) {
      return;
    }

    for (const auto& muon : muons) {
      if (muon.trackType() != mTrackType) {
        continue;
      }
      if (mApplyMcMask && (muon.mcMask() != 0)) {
        continue;
      }
      const auto eta(muon.eta());
      if ((eta >= -2.5) || (eta < -3.6)) {
        continue;
      }

      fillHistograms(muon);
    } // loop over muons

    return;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskMuonSourceMC>(cfgc),
  };
}
