// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Configurable.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetMatching.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MCTruthContainer.h"

#include <cmath>

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(NConst, nConst, int);
DECLARE_SOA_COLUMN(Girth, girth, float);
DECLARE_SOA_COLUMN(PTD, pTD, float);
DECLARE_SOA_COLUMN(MatchDeltaR, matchDeltaR, float);
DECLARE_SOA_COLUMN(PtResponse, ptResponse, float);
DECLARE_SOA_COLUMN(QGLabel, qgLabel, int);

DECLARE_SOA_TABLE(QGJetTable, "AOD", "QGJET",
                  JetPt,
                  JetEta,
                  JetPhi,
                  NConst,
                  Girth,
                  PTD,
                  MatchDeltaR,
                  PtResponse,
                  QGLabel);
}

//------------------------------------------------
// helper functions
//------------------------------------------------
float deltaPhi(float phi1, float phi2)
{
  return std::remainder(phi1 - phi2, 2.f * static_cast<float>(M_PI));
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = deltaPhi(phi1, phi2);
  return std::sqrt(deta * deta + dphi * dphi);
}

//------------------------------------------------
// find initiating parton by ancestry
//------------------------------------------------
int getInitiatingParton(auto const& particle,
                        aod::McParticles const& mcParticles)
{
  auto p = particle;
  int pdg = p.pdgCode();

  while (p.has_mothers()) {
    auto mothers = p.mothers_as<aod::McParticles>();
    if (mothers.size() == 0) {
      break;
    }

    auto mom = mothers.iteratorAt(0);
    int mpdg = mom.pdgCode();

    // stop at quark or gluon
    if (std::abs(mpdg) == 21 || (std::abs(mpdg) >= 1 && std::abs(mpdg) <= 6)) {
      pdg = mpdg;
    }

    p = mom;
  }

  return pdg;
}

//------------------------------------------------
// main task
//------------------------------------------------
struct QGTreeCreator {

  Configurable<float> jetPtMin{"jetPtMin",10.f};
  Configurable<float> maxMatchDeltaR{"maxMatchDeltaR",0.3f};

  Produces<aod::QGJetTable> qgjets;

  void process(aod::ChargedMCDetectorLevelJets const& recoJets,
               aod::ChargedMCParticleLevelJets const& truthJets,
               aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets const& matches,
               aod::McParticles const& mcParticles)
  {
    for (auto const& jet : recoJets) {

      if (jet.pt() < jetPtMin)
        continue;

      //----------------------------------
      // compute jet observables
      //----------------------------------
      int nconst = 0;
      float sumPt = 0;
      float sumPt2 = 0;
      float sumPtDr = 0;

      for (auto const& c : jet.tracks_as<aod::ChargedMCDetectorLevelJetConstituent>()) {
        float pt = c.pt();
        float dr = deltaR(c.eta(), c.phi(), jet.eta(), jet.phi());

        nconst++;
        sumPt += pt;
        sumPt2 += pt*pt;
        sumPtDr += pt*dr;
      }

      float girth = sumPt>0 ? sumPtDr/sumPt : -1;
      float ptd   = sumPt>0 ? std::sqrt(sumPt2)/sumPt : -1;

      //----------------------------------
      // matching block
      //----------------------------------
      float matchDr = -1;
      float ptResp = -1;
      int qg = -1;

      for (auto const& match : matches) {

        if (match.chargedMCDetectorLevelJetId() != jet.globalIndex())
          continue;

        auto truthJet = truthJets.iteratorAt(
            match.chargedMCParticleLevelJetId());

        matchDr = deltaR(jet.eta(), jet.phi(),
                         truthJet.eta(), truthJet.phi());

        if (matchDr > maxMatchDeltaR)
          continue;

        ptResp = jet.pt() / truthJet.pt();

        //----------------------------------
        // find initiating parton
        //----------------------------------
        float maxPt = -1;
        int pdg = 0;

        for (auto const& tc :
             truthJet.tracks_as<aod::ChargedMCParticleLevelJetConstituent>())
        {
          if (!tc.has_mcParticle())
            continue;

          auto mc = tc.mcParticle();

          if (tc.pt() > maxPt) {
            maxPt = tc.pt();
            pdg = getInitiatingParton(mc, mcParticles);
          }
        }

        //----------------------------------
        // assign q/g label
        //----------------------------------
        if (std::abs(pdg) == 21)
          qg = 1; // gluon
        else if (std::abs(pdg) >= 1 && std::abs(pdg) <= 6)
          qg = 0; // quark

        break;
      }

      //----------------------------------
      // store
      //----------------------------------
      qgjets(jet.pt(),
             jet.eta(),
             jet.phi(),
             nconst,
             girth,
             ptd,
             matchDr,
             ptResp,
             qg);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QGTreeCreator>(cfgc, TaskName{"qg-tree-creator"})};
}
