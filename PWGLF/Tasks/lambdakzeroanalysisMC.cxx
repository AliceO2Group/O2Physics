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
// Adaptation of V0 analysis task for run on MC data
// ========================
//
// This code loops over a V0Data table and produces some
// standard analysis output. It requires either
// the lambdakzerofinder or the lambdakzeroproducer tasks
// to have been executed in the workflow (before).
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    aimeric.landou@cern.ch (MC adaptation)
//    david.dobrigkeit.chinellato@cern.ch (original lambdakzeroanalysis task)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/MC.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExt = soa::Join<aod::FullTracks, aod::TracksExtended>;
using MyTracks = soa::Join<FullTracksExt, aod::McTrackLabels>;

struct lambdakzeroQA {
  //Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassLambda", "hMassLambda", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},

      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{2000, -10.0f, 10.0f}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{2000, -10.0f, 10.0f}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
    },
  };

  void process(aod::Collision const& collision, aod::V0Datas const& fullV0s)
  {
    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());
      registry.fill(HIST("hMassLambda"), v0.mLambda());
      registry.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());

      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());
    }
  }
};

struct lambdakzeroanalysisMC {

  HistogramRegistry registry{
    "registry",
    {
      {"h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 0.450f, 0.550f}}}},
      {"h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"h3dMassK0Short_MCportion", "h3dMassK0Short_MCportion", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 0.450f, 0.550f}}}},
      {"h3dMassLambda_MCportion", "h3dMassLambda_MCportion", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"h3dMassAntiLambda_MCportion", "h3dMassAntiLambda_MCportion", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"MCmomID_Lambda", "MCmomID_Lambda", {HistType::kTH1I, {{4000000, 0, 4000000}}}},
      {"MCmomID_AntiLambda", "MCmomID_AntiLambda", {HistType::kTH1I, {{4000000, 0, 4000000}}}},
      {"MCmomID_K0Short", "MCmomID_K0Short", {HistType::kTH1I, {{4000000, 0, 4000000}}}},
      {"V0loopFiltersCounts", "V0loopFiltersCounts", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},
    },
  };
  //_MCportion: taking the h3dMassXXX histogram and removing all v0s for which the postrack and negtrack mothers' don't match or are not K0S

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};
  ConfigurableAxis massK0Shortbinning{"K0S-mass-binning", {200, 0.450f, 0.550f}, ""};
  ConfigurableAxis massLambdabinning{"Lambda-mass-binning", {200, 1.015f, 1.215f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisK0Short = {massK0Shortbinning, "Inv. Mass (GeV)"};
    AxisSpec massAxisLambda = {massLambdabinning, "Inv. Mass (GeV)"};

    registry.add("h3dMassK0ShortDca", "h3dMassK0ShortDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisK0Short}});
    registry.add("h3dMassLambdaDca", "h3dMassLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
    registry.add("h3dMassAntiLambdaDca", "h3dMassAntiLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
    registry.add("h3dMassK0ShortDca_MCportion", "h3dMassK0ShortDca_MCportion", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisK0Short}});
    registry.add("h3dMassLambdaDca_MCportion", "h3dMassLambdaDca_MCportion", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
    registry.add("h3dMassAntiLambdaDca_MCportion", "h3dMassAntiLambdaDca_MCportion", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});

    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(1, "V0 Candidates");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(2, "V0Radius and CosPA");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(4, "Lambda Rapidity");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(5, "Lambda lifetime cut");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(7, "K0S Rapidity");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(8, "K0S lifetime cut");
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};
  Configurable<int> checksel7{"checksel7", 0, "checksel7"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  //   void process(aod::Collision const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    // if (!collision.alias()[kINT7]) { //should not use if looking at MC
    //   return;
    // }
    if (checksel7 == 1 && !collision.sel7()) {
      return;
    }

    for (auto& v0 : fullV0s) {
      //   FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("V0loopFiltersCounts"), 0.5);
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("V0loopFiltersCounts"), 1.5);

        auto mcnegtrack = v0.negTrack_as<MyTracks>().mcParticle();
        auto mcpostrack = v0.posTrack_as<MyTracks>().mcParticle();
        auto particleMotherOfNeg = mcnegtrack.mother0_as<aod::McParticles>();
        bool MomIsPrimary = MC::isPhysicalPrimary(particleMotherOfNeg);

        if (TMath::Abs(v0.yLambda()) < rapidity) {
          registry.fill(HIST("V0loopFiltersCounts"), 3.5);
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0) < lifetimecut->get("lifetimecutLambda")) {
            registry.fill(HIST("V0loopFiltersCounts"), 4.5);
            registry.fill(HIST("h3dMassLambda"), collision.centV0M(), v0.pt(), v0.mLambda());
            registry.fill(HIST("h3dMassAntiLambda"), collision.centV0M(), v0.pt(), v0.mAntiLambda());

            if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == 3122) {
              registry.fill(HIST("h3dMassLambda_MCportion"), collision.centV0M(), v0.pt(), v0.mLambda());
              registry.fill(HIST("MCmomID_Lambda"), mcnegtrack.mother0Id());
            }
            if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == -3122) {
              registry.fill(HIST("h3dMassAntiLambda_MCportion"), collision.centV0M(), v0.pt(), v0.mAntiLambda());
              registry.fill(HIST("MCmomID_AntiLambda"), mcnegtrack.mother0Id());
            }
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mLambda());
              registry.fill(HIST("h3dMassAntiLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiLambda());

              if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == 3122) {
                registry.fill(HIST("h3dMassLambdaDca_MCportion"), v0.dcaV0daughters(), v0.pt(), v0.mLambda());
              }
              if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == -3122) {
                registry.fill(HIST("h3dMassAntiLambdaDca_MCportion"), v0.dcaV0daughters(), v0.pt(), v0.mAntiLambda());
              }
            }
          }
        }
        if (TMath::Abs(v0.yK0Short()) < rapidity) {
          registry.fill(HIST("V0loopFiltersCounts"), 6.5);
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short) < lifetimecut->get("lifetimecutK0S")) {
            registry.fill(HIST("V0loopFiltersCounts"), 7.5);

            registry.fill(HIST("h3dMassK0Short"), collision.centV0M(), v0.pt(), v0.mK0Short());

            if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == 310) {
              registry.fill(HIST("h3dMassK0Short_MCportion"), collision.centV0M(), v0.pt(), v0.mK0Short());
              registry.fill(HIST("MCmomID_K0Short"), mcnegtrack.mother0Id());
            }
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassK0ShortDca"), v0.dcaV0daughters(), v0.pt(), v0.mK0Short());
              if (MomIsPrimary && mcnegtrack.mother0Id() == mcpostrack.mother0Id() && particleMotherOfNeg.pdgCode() == 310) {
                registry.fill(HIST("h3dMassK0ShortDca_MCportion"), v0.dcaV0daughters(), v0.pt(), v0.mK0Short());
              }
            }
          }
        }
      }
    }
  }
};

struct lambdakzeroParticleCountMC {
  //Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hK0ShortCount", "hK0ShortCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hLambdaCount", "hLambdaCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hAntiLambdaCount", "hAntiLambdaCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hK0ShortCount_PtDiff", "hK0ShortCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hLambdaCount_PtDiff", "hLambdaCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hAntiLambdaCount_PtDiff", "hAntiLambdaCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
    },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hK0ShortCount"))->GetXaxis()->SetBinLabel(1, "K0S mothers");
    registry.get<TH1>(HIST("hK0ShortCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
    registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(1, "Lambda mothers");
    registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
    registry.get<TH1>(HIST("hAntiLambdaCount"))->GetXaxis()->SetBinLabel(1, "AntiLambda mothers");
    registry.get<TH1>(HIST("hAntiLambdaCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
  }

  Configurable<float> rapidityMCcut{"rapidityMCcut", 0.5, "rapidityMCcut"};

  // void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  void process(aod::McParticles const& mcParticles)
  {
    for (auto& mcparticle : mcParticles) {
      if (TMath::Abs(mcparticle.y()) < rapidityMCcut) {
        if (MC::isPhysicalPrimary(mcparticle)) {
          if (mcparticle.pdgCode() == 310) {
            registry.fill(HIST("hK0ShortCount"), 0.5);
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -211) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 211)) {
              registry.fill(HIST("hK0ShortCount"), 1.5);
              registry.fill(HIST("hK0ShortCount_PtDiff"), mcparticle.pt());
            }
          }
          if (mcparticle.pdgCode() == 3122) {
            registry.fill(HIST("hLambdaCount"), 0.5);
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 2212) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 2212 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -211)) {
              registry.fill(HIST("hLambdaCount"), 1.5);
              registry.fill(HIST("hLambdaCount_PtDiff"), mcparticle.pt());
            }
          }
          if (mcparticle.pdgCode() == -3122) {
            registry.fill(HIST("hAntiLambdaCount"), 0.5);
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -2212) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -2212 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 211)) {
              registry.fill(HIST("hAntiLambdaCount"), 1.5);
              registry.fill(HIST("hAntiLambdaCount_PtDiff"), mcparticle.pt());
            }
          }
        }
      }
    }
  }
};

struct V0daughtersTrackingEfficiency {

  HistogramRegistry registry{
    "registry",
    {
      {"hV0DaughterMcParticleIDs", "hV0DaughterMcParticleIDs", {HistType::kTH1I, {{1000000, 0, 1000000}}}},
      {"hAssociatedMcParticleIDs", "hAssociatedMcParticleIDs", {HistType::kTH1I, {{1000000, 0, 1000000}}}},
    },
  };

  void process(aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    for (auto& mcparticle : mcParticles) {
      if (MC::isPhysicalPrimary(mcparticle)) {
        if (TMath::Abs(mcparticle.y()) < 0.5) {
          if (mcparticle.pdgCode() == 310) {
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -211) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 211)) {
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter0Id());
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter1Id());
            }
          }
          if (mcparticle.pdgCode() == 3122) {
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 2212) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 2212 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -211)) {
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter0Id());
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter1Id());
            }
          }
          if (mcparticle.pdgCode() == -3122) {
            if ((mcparticle.daughter0_as<aod::McParticles>().pdgCode() == 211 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == -2212) || (mcparticle.daughter0_as<aod::McParticles>().pdgCode() == -2212 && mcparticle.daughter1_as<aod::McParticles>().pdgCode() == 211)) {
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter0Id());
              registry.fill(HIST("hV0DaughterMcParticleIDs"), mcparticle.daughter1Id());
            }
          }
        }
      }
    }
    for (auto& track : tracks) {
      auto AssociatedMcParticle = track.mcParticle();
      if (TMath::Abs(AssociatedMcParticle.y()) < 0.5) {
        registry.fill(HIST("hAssociatedMcParticleIDs"), AssociatedMcParticle.globalIndex());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroanalysisMC>(cfgc, TaskName{"lf-lambdakzeroanalysisMC"}),
    adaptAnalysisTask<lambdakzeroQA>(cfgc, TaskName{"lf-lambdakzeroQA"}),
    adaptAnalysisTask<lambdakzeroParticleCountMC>(cfgc, TaskName{"lf-lambdakzeroParticleCountMC"}),
    adaptAnalysisTask<V0daughtersTrackingEfficiency>(cfgc, TaskName{"lf-V0daughtersTrackingEfficiency"})};
}
