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

#include <Math/Vector4D.h>
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
using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
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
  Configurable<int> trackType{"trackType", 3, "Muon track type, validated values are 0, 1, 2, 3 and 4"};
  Configurable<int> charge{"charge", 0, "Muon track charge, validated values are 0, 1 and -1, 0 represents both 1 and -1"};
  Configurable<bool> pairSource{"pairSource", true, "check also the source of like-sign muon pairs"};

  double pDcaMax = 594.0;  // p*DCA maximum value for small Rab
  double pDcaMax2 = 324.0; // p*DCA maximum value for large Rabs
  double rAbsMid = 26.5;   // R at absorber end minimum value
  double rAbsMax = 89.5;   // R at absorber end maximum value
  double rAbsMin = 17.6;   // R at absorber end maximum value
  double etaLow = -4.0;    // low edge of eta acceptance
  double etaUp = -2.5;     // up edge of eta acceptance
  double edgeZ = 10.0;     // edge of event position Z
  double ptLow = 1.0;      // low edge of pT for muon pairs

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
    AxisSpec const axisMass{200, 0., 20., "Inv.Mass (GeV/#it{c}^2)"};

    HistogramConfigSpec const h1ColNumber{HistType::kTH1F, {axisColNumber}};
    HistogramConfigSpec const h1Pt{HistType::kTH1F, {axisPt}};
    HistogramConfigSpec h1Mass{HistType::kTH1F, {axisMass}};
    HistogramConfigSpec const h2PtDCA{HistType::kTH2F, {axisPt, axisDCA}};
    HistogramConfigSpec const h2PtChi2{HistType::kTH2F, {axisPt, axisChi2}};
    HistogramConfigSpec const h2PtDeltaPt{HistType::kTH2F, {axisPt, axisDeltaPt}};

    registry.add("h1ColNumber", "", h1ColNumber);
    registry.add("h1MuBeforeCuts", "", h1Pt);
    registry.add("h1MuonMass", "", h1Mass);
    registry.add("h1BeautyMass", "", h1Mass);
    registry.add("h1OtherMass", "", h1Mass);
    registry.add("h1MuonMassGen", "", h1Mass);
    registry.add("h1BeautyMassGen", "", h1Mass);
    registry.add("h1OtherMassGen", "", h1Mass);
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
      if (pdgAbs < 10 || pdgAbs == 21) {
        break; // Quark and gluon
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
      if ((pdgRem % 100 == 1 || pdgRem % 100 == 3) && pdgRem > 1000) { // diquarks
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

  int traceAncestor(const McMuons::iterator& muon, aod::McParticles const& mctracks)
  {
    int mcNum = 0;
    if (!muon.has_mcParticle()) {
      return 0;
    }
    auto mcPart(muon.mcParticle());
    if (std::abs(mcPart.pdgCode()) != kMuonMinus) {
      return 0;
    }
    while (mcPart.has_mothers()) { // the first hadron after hadronization
      auto mother = mcPart.mothers_first_as<aod::McParticles>();
      if (std::abs(mother.getGenStatusCode()) < 80) {
        break;
      }
      mcPart = mother;
    }
    int flv = mcPart.pdgCode() / std::pow(10, static_cast<int>(std::log10(std::abs(mcPart.pdgCode()))));
    if (abs(flv) == 5 && mcPart.pdgCode() < 1000)
      flv = -flv;
    for (int i = (mcPart.mothers_first_as<aod::McParticles>()).globalIndex(); i <= (mcPart.mothers_last_as<aod::McParticles>()).globalIndex(); i++) { // loop over the lund string
      for (auto mctrack : mctracks) {
        if (mctrack.globalIndex() != i) {
          continue;
        }
        if ((mctrack.pdgCode() != flv) && (abs(mctrack.pdgCode()) < abs(flv) * 1000)) {
          continue;
        }
        while (mctrack.has_mothers()) {
          int motherflv = (mctrack.mothers_first_as<aod::McParticles>()).pdgCode() / std::pow(10, static_cast<int>(std::log10(abs((mctrack.mothers_first_as<aod::McParticles>()).pdgCode())))); // find the mother with same flavor
          auto mother = (abs(motherflv) == abs(flv)) ? (mctrack.mothers_first_as<aod::McParticles>()) : (mctrack.mothers_last_as<aod::McParticles>());
          if ((mother.pdgCode() != mctrack.pdgCode()) && (abs(mctrack.pdgCode()) < 10)) { // both mother is not the the quark with same flavor
            mcNum = mctrack.globalIndex();
            return mcNum;
          }
          mctrack = mother;
        }
      }
    }
    return 0;
  }
  bool Corr(const McMuons::iterator& muon1, const McMuons::iterator& muon2, aod::McParticles const& mcParts)
  {

    int moth11(0), moth12(0), moth21(1), moth22(1);
    int anc1 = traceAncestor(muon1, mcParts);
    int anc2 = traceAncestor(muon2, mcParts);
    if (anc1 == 0 || anc2 == 0) {
      return false;
    }
    for (auto mcPart : mcParts) {
      if (mcPart.globalIndex() == anc1) {
        moth11 = (mcPart.mothers_first_as<aod::McParticles>()).globalIndex();
        moth12 = (mcPart.mothers_last_as<aod::McParticles>()).globalIndex();
      }
      if (mcPart.globalIndex() == anc2) {
        moth21 = (mcPart.mothers_first_as<aod::McParticles>()).globalIndex();
        moth22 = (mcPart.mothers_last_as<aod::McParticles>()).globalIndex();
      }
    }
    if ((moth11 == moth21) && (moth12 == moth22)) {
      return true;
    }
    return false; // uncorrelated
  }
  void fillPairs(const McMuons::iterator& muon, const McMuons::iterator& muon2, aod::McParticles const& mcParts)
  {
    if (trackType != 3) {
      return;
    }
    float mm = o2::constants::physics::MassMuon;

    const auto mask1(getMask(muon));
    const auto mask2(getMask(muon2));

    ROOT::Math::PtEtaPhiMVector mu1Vec(muon.pt(), muon.eta(), muon.phi(), mm);
    ROOT::Math::PtEtaPhiMVector mu2Vec(muon2.pt(), muon2.eta(), muon2.phi(), mm);
    ROOT::Math::PtEtaPhiMVector dimuVec = mu1Vec + mu2Vec;
    auto InvM = dimuVec.M();

    if (!muon.has_mcParticle() || !muon2.has_mcParticle()) {
      return;
    }
    auto mcPart1(muon.mcParticle());
    auto mcPart2(muon2.mcParticle());

    ROOT::Math::PtEtaPhiMVector mu1VecGen(mcPart1.pt(), mcPart1.eta(), mcPart1.phi(), mm);
    ROOT::Math::PtEtaPhiMVector mu2VecGen(mcPart2.pt(), mcPart2.eta(), mcPart2.phi(), mm);
    ROOT::Math::PtEtaPhiMVector dimuVecGen = mu1VecGen + mu2VecGen;
    auto InvMGen = dimuVecGen.M();

    if (isMuon(mask1) && isMuon(mask2)) {
      registry.fill(HIST("h1MuonMass"), InvM);
      registry.fill(HIST("h1MuonMassGen"), InvMGen);
    }
    if (Corr(muon, muon2, mcParts) && isBeautyMu(mask1) && isBeautyMu(mask2)) {
      registry.fill(HIST("h1BeautyMass"), InvM);
      registry.fill(HIST("h1BeautyMassGen"), InvMGen);
    } else {
      registry.fill(HIST("h1OtherMass"), InvM);
      registry.fill(HIST("h1OtherMassGen"), InvMGen);
    }
  }

  void process(MyCollisions::iterator const& collision,
               McMuons const& muons,
               aod::McParticles const& mcParts)
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
      if ((rAbs >= rAbsMid) || (rAbs < rAbsMin)) {
        if (pDca >= pDcaMax || pDca < 0) {
          continue;
        }
      }
      if ((rAbs >= rAbsMax) || (rAbs < rAbsMid)) {
        if (pDca >= pDcaMax2 || pDca < 0) {
          continue;
        }
      }
      if ((muon.chi2() >= 1e6) || (muon.chi2() < 0)) {
        continue;
      }
      if (charge != 0 && muon.sign() != charge) {
        continue;
      }
      fillHistograms(muon);
      if (pairSource) {
        if (muon.pt() < ptLow) {
          continue;
        }
        for (const auto& muon2 : muons) {
          if (muon2.sign() != muon.sign()) {
            continue;
          }
          if (muon2.globalIndex() <= muon.globalIndex()) {
            continue;
          }
          // muon selections
          if (muon2.trackType() != trackType) {
            continue;
          }
          if (muon2.pt() < ptLow) {
            continue;
          }
          const auto eta2(muon2.eta()), pDca2(muon2.pDca()), rAbs2(muon2.rAtAbsorberEnd());
          if ((eta2 >= etaUp) || (eta2 < etaLow)) {
            continue;
          }
          if ((rAbs2 >= rAbsMid) || (rAbs2 < rAbsMin)) {
            if (pDca2 >= pDcaMax || pDca2 < 0) {
              continue;
            }
          }
          if ((rAbs2 >= rAbsMax) || (rAbs2 < rAbsMid)) {
            if (pDca2 >= pDcaMax2 || pDca2 < 0) {
              continue;
            }
          }

          if ((muon2.chi2() >= 1e6) || (muon2.chi2() < 0)) {
            continue;
          }
          if ((muon2.chi2MatchMCHMID() >= 1e6) || (muon2.chi2MatchMCHMID() < 0)) {
            continue;
          }
          fillPairs(muon, muon2, mcParts);
        }
      }
    } // loop over muons
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleMuonSource>(cfgc),
  };
}
