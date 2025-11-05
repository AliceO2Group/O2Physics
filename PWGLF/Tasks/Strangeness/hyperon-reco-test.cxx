// Copyright 2020-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TMath.h"
#include <TDatabasePDG.h>
#include <TH2F.h>
#include <TPDGCode.h>

#include <cmath>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct myLambda {
  HistogramRegistry registry{
    "registry",
    {{"hReduction", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hReduction_mc", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hmyLambda_before",
      "",
      {HistType::kTH2F, {{250, 1.0, 1.25, "Lambda_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcLambda_before",
      "",
      {HistType::kTH2F, {{250, 1.0, 1.25, "Lambda_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyLambda_after",
      "",
      {HistType::kTH2F, {{250, 1.0, 1.25, "Lambda_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcLambda_after",
      "",
      {HistType::kTH2F, {{250, 1.0, 1.25, "Lambda_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyLambda_after_b", "", {HistType::kTH1F, {{250, 1.0, 1.25}}}},
     {"hmcLambda_after_b", "", {HistType::kTH1F, {{250, 1.0, 1.25}}}}}};

  // track cuts
  Configurable<float> rapidity{"rapidity", 0.80, "V0 rapidity"};
  Configurable<float> tpcNcl{"tpcNcl", 70, "tpcNcl"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> tpcsigma{"tpcsigma", 5, "tpcsigma"};

  // V0 cuts
  Configurable<float> minpt{"minpt", 0.30, "minpt"};
  Configurable<float> v0radius{"v0radius", 0.2, "v0radius"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA Pos To PV"};
  Configurable<float> v0cospa{"v0cospa", 0.99, "V0 CosPA"};
  Configurable<float> removeKs{"removeKs", 0.025, "removeKs"};

  Filter filtertrack = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& /*collision*/,
               soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s,
               DauTracks const& /*tracks*/,
               aod::McParticles const& /*mcParticles*/)
  {
    for (auto& v0 : V0s) {
      registry.fill(HIST("hReduction"), 0.5);
      registry.fill(HIST("hmyLambda_before"), v0.mLambda(), v0.pt());

      if (v0.has_mcParticle()) {
        auto v0mcparticle = v0.mcParticle();
        if (v0mcparticle.pdgCode() == 3122) {
          registry.fill(HIST("hReduction_mc"), 0.5);
          registry.fill(HIST("hmcLambda_before"), v0.mLambda(), v0.pt());
        }
      }

      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();

      if ( // track cuts
        TMath::Abs(posdau.eta()) < rapidity && TMath::Abs(negdau.eta()) < rapidity && posdau.tpcNClsCrossedRows() > tpcNcl && negdau.tpcNClsCrossedRows() > tpcNcl && TMath::Abs(posdau.tpcNSigmaPr()) < tpcsigma && TMath::Abs(negdau.tpcNSigmaPi()) < tpcsigma &&

        // V0 cuts
        v0.pt() > minpt && v0.v0radius() > v0radius && v0.v0cosPA() > v0cospa && TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) > removeKs) {
        registry.fill(HIST("hReduction"), 1.5);
        registry.fill(HIST("hmyLambda_after"), v0.mLambda(), v0.pt());
        registry.fill(HIST("hmyLambda_after_b"), v0.mLambda());

        if (v0.has_mcParticle()) {
          auto v0mcparticle = v0.mcParticle();
          if (v0mcparticle.pdgCode() == 3122) {
            registry.fill(HIST("hReduction_mc"), 1.5);
            registry.fill(HIST("hmcLambda_after"), v0.mLambda(), v0.pt());
            registry.fill(HIST("hmcLambda_after_b"), v0.mLambda());
          }
        }
      }
    }
  }
};

struct myXi {
  HistogramRegistry registry{
    "registry",
    {{"hReduction", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hReduction_mc", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hmyXi_before",
      "",
      {HistType::kTH2F, {{500, 1.2, 1.7, "Xi_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcXi_before",
      "",
      {HistType::kTH2F, {{500, 1.2, 1.7, "Xi_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyLambda", "", {HistType::kTH1F, {{250, 1.0, 1.25}}}},
     {"hmyXi_after", "", {HistType::kTH2F, {{500, 1.2, 1.7, "Xi_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcXi_after", "", {HistType::kTH2F, {{500, 1.2, 1.7, "Xi_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyXi_after_b", "", {HistType::kTH1F, {{500, 1.2, 1.7}}}},
     {"hmcXi_after_b", "", {HistType::kTH1F, {{500, 1.2, 1.7}}}}}};

  // Dau & bach track cuts
  Configurable<float> rapidity{"rapidity", 0.80, "V0 rapidity"};
  Configurable<float> mincrossedrow{"mincrossedrow", 70, "mincrossedrow"};
  Configurable<float> minpt{"minpt", 0.30, "minpt"};
  Configurable<float> dcatopv{"dcatopv", 0.05, "DCA To PV"};
  Configurable<float> tpcsigma{"tpcsigma", 4.0, "tpcsigma"};

  // V0 cuts
  Configurable<float> v0radius{"v0radius", 1.4, "v0radius"};
  Configurable<float> dcav0dau{"dcav0dau", 1.6, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 0.07, "dcav0topv"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "v0CPA"};
  Configurable<float> v0masswindow{"v0masswindow", 0.006, "v0masswindow"};

  // casc cuts
  Configurable<float> cascradius{"cascradius", 0.8, "cascradius"};
  Configurable<float> dcacascdau{"dcascascdau", 1.6, "DCA V0 to bach"};
  Configurable<float> casccospa{"casccospa", 0.98, "casc CosPA"};
  Configurable<float> removeOmega{"removeXi", 0.005, "removeXi"};

  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  void process(aod::Collision const& collision,
               soa::Join<aod::CascDataExt, aod::McCascLabels> const& Cascades,
               DauTracks const& /*tracks*/,
               aod::McParticles const& /*mcParticles*/)
  {
    for (auto& casc : Cascades) {
      registry.fill(HIST("hReduction"), 0.5);
      registry.fill(HIST("hmyXi_before"), casc.mXi(), casc.pt());

      if (casc.has_mcParticle()) {
        auto cascmcparticle = casc.mcParticle();
        if (cascmcparticle.pdgCode() == 3312) {
          registry.fill(HIST("hReduction_mc"), 0.5);
          registry.fill(HIST("hmcXi_before"), casc.mXi(), casc.pt());
        }
      }

      auto posdau = casc.posTrack_as<DauTracks>();
      auto negdau = casc.negTrack_as<DauTracks>();
      auto bachelor = casc.bachelor_as<DauTracks>();

      if ( // Dau & bach track cuts
        TMath::Abs(posdau.eta()) < rapidity && TMath::Abs(negdau.eta()) < rapidity && TMath::Abs(bachelor.eta()) < rapidity && posdau.tpcNClsCrossedRows() > mincrossedrow && negdau.tpcNClsCrossedRows() > mincrossedrow && bachelor.tpcNClsCrossedRows() > mincrossedrow && posdau.pt() > minpt && negdau.pt() > minpt && bachelor.pt() > minpt && TMath::Abs(casc.dcapostopv()) > dcatopv && TMath::Abs(casc.dcanegtopv()) > dcatopv && TMath::Abs(casc.dcabachtopv()) > dcatopv && TMath::Abs(posdau.tpcNSigmaPr()) < tpcsigma && TMath::Abs(negdau.tpcNSigmaPi()) < tpcsigma && TMath::Abs(bachelor.tpcNSigmaPi()) < tpcsigma && bachelor.sign() < 0 &&

        // V0 cuts
        casc.v0radius() > v0radius && casc.dcaV0daughters() < dcav0dau && casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv && casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("hmyLambda"), casc.mLambda());

        if (TMath::Abs(casc.mLambda() - o2::constants::physics::MassLambda0) < v0masswindow &&

            // Cascade cut
            casc.cascradius() > cascradius && casc.dcacascdaughters() < dcacascdau && casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa && TMath::Abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > removeOmega && casc.sign() < 0) {
          registry.fill(HIST("hReduction"), 1.5);
          registry.fill(HIST("hmyXi_after"), casc.mXi(), casc.pt());
          registry.fill(HIST("hmyXi_after_b"), casc.mXi());

          if (casc.has_mcParticle()) {
            auto cascmcparticle = casc.mcParticle();
            if (cascmcparticle.pdgCode() == 3312) {
              registry.fill(HIST("hReduction_mc"), 1.5);
              registry.fill(HIST("hmcXi_after"), casc.mXi(), casc.pt());
              registry.fill(HIST("hmcXi_after_b"), casc.mXi());
            }
          }
        }
      }
    }
  }
};

struct myOmega {
  HistogramRegistry registry{
    "registry",
    {{"hReduction", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hReduction_mc", "", {HistType::kTH1F, {{2, 0, 2}}}},
     {"hmyOmega_before",
      "",
      {HistType::kTH2F, {{500, 1.5, 2.0, "Omega_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcOmega_before",
      "",
      {HistType::kTH2F, {{500, 1.5, 2.0, "Omega_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyLambda", "", {HistType::kTH1F, {{250, 1.0, 1.25}}}},
     {"hmyOmega_after",
      "",
      {HistType::kTH2F, {{500, 1.5, 2.0, "Omega_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmcOmega_after",
      "",
      {HistType::kTH2F, {{500, 1.5, 2.0, "Omega_mass"}, {10000, 0, 10, "pt"}}}},
     {"hmyOmega_after_b", "", {HistType::kTH1F, {{500, 1.5, 2.0}}}},
     {"hmcOmega_after_b", "", {HistType::kTH1F, {{500, 1.5, 2.0}}}}}};

  // Dau & bach track cuts
  Configurable<float> rapidity{"rapidity", 0.80, "V0 rapidity"};
  Configurable<float> mincrossedrow{"mincrossedrow", 70, "mincrossedrow"};
  Configurable<float> minpt{"minpt", 0.30, "minpt"};
  Configurable<float> dcatopv{"dcatopv", 0.05, "DCA To PV"};
  Configurable<float> tpcsigma{"tpcsigma", 4.0, "tpcsigma"};

  // V0 cuts
  Configurable<float> v0radius{"v0radius", 1.0, "v0radius"};
  Configurable<float> dcav0dau{"dcav0dau", 1.2, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 0.06, "dcav0topv"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "v0CPA"};
  Configurable<float> v0masswindow{"v0masswindow", 0.006, "v0masswindow"};

  // casc cuts
  Configurable<float> cascradius{"cascradius", 0.2, "cascradius"};
  Configurable<float> dcacascdau{"dcascascdau", 0.8, "DCA V0 to bach"};
  Configurable<float> casccospa{"casccospa", 0.995, "casc CosPA"};
  Configurable<float> removeXi{"removeXi", 0.005, "removeXi"};

  using DauTracks =
    soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

  void process(aod::Collision const& collision,
               soa::Join<aod::CascDataExt, aod::McCascLabels> const& Cascades,
               DauTracks const& /*tracks*/,
               aod::McParticles const& /*mcParticles*/)
  {
    for (auto& casc : Cascades) {
      registry.fill(HIST("hReduction"), 0.5);
      registry.fill(HIST("hmyOmega_before"), casc.mOmega(), casc.pt());

      if (casc.has_mcParticle()) {
        auto cascmcparticle = casc.mcParticle();
        if (cascmcparticle.pdgCode() == 3334) {
          registry.fill(HIST("hReduction_mc"), 0.5);
          registry.fill(HIST("hmcOmega_before"), casc.mOmega(), casc.pt());
        }
      }

      auto posdau = casc.posTrack_as<DauTracks>();
      auto negdau = casc.negTrack_as<DauTracks>();
      auto bachelor = casc.bachelor_as<DauTracks>();

      if ( // Dau & bach track cuts
        TMath::Abs(posdau.eta()) < rapidity && TMath::Abs(negdau.eta()) < rapidity && TMath::Abs(bachelor.eta()) < rapidity && posdau.tpcNClsCrossedRows() > mincrossedrow && negdau.tpcNClsCrossedRows() > mincrossedrow && bachelor.tpcNClsCrossedRows() > mincrossedrow && posdau.pt() > minpt && negdau.pt() > minpt && bachelor.pt() > minpt && TMath::Abs(casc.dcapostopv()) > dcatopv && TMath::Abs(casc.dcanegtopv()) > dcatopv && TMath::Abs(casc.dcabachtopv()) > dcatopv && TMath::Abs(posdau.tpcNSigmaPr()) < tpcsigma && TMath::Abs(negdau.tpcNSigmaPi()) < tpcsigma && TMath::Abs(bachelor.tpcNSigmaKa()) < tpcsigma && bachelor.sign() < 0 &&

        // V0 cuts
        casc.v0radius() > v0radius && casc.dcaV0daughters() < dcav0dau && casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv && casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("hmyLambda"), casc.mLambda());

        if (TMath::Abs(casc.mLambda() - o2::constants::physics::MassLambda0) < v0masswindow &&

            // Cascade cut
            casc.cascradius() > cascradius && casc.dcacascdaughters() < dcacascdau && casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa && TMath::Abs(casc.mXi() - o2::constants::physics::MassXiMinus) > removeXi && casc.sign() < 0) {
          registry.fill(HIST("hReduction"), 1.5);
          registry.fill(HIST("hmyOmega_after"), casc.mOmega(), casc.pt());
          registry.fill(HIST("hmyOmega_after_b"), casc.mOmega());

          if (casc.has_mcParticle()) {
            auto cascmcparticle = casc.mcParticle();
            if (cascmcparticle.pdgCode() == 3334) {
              registry.fill(HIST("hReduction_mc"), 1.5);
              registry.fill(HIST("hmcOmega_after"), casc.mOmega(), casc.pt());
              registry.fill(HIST("hmcOmega_after_b"), casc.mOmega());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<myLambda>(cfgc),
                      adaptAnalysisTask<myXi>(cfgc),
                      adaptAnalysisTask<myOmega>(cfgc)};
}
