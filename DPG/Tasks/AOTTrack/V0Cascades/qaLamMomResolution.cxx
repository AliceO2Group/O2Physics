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

///
/// \file   qaLamMomResolution.cxx
/// \author Carolina Reetz c.reetz@cern.ch
/// \brief  QA task to study momentum resolution of Lambda daughter tracks
/// \TODO:  differentiate between pos and neg eta

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using V0DatasLabeled = soa::Join<aod::V0Datas, aod::V0Covs, aod::McV0Labels>;

struct qaLamMomResolution {
  ConfigurableAxis mBins{"mBins", {200, 1.0f, 1.2f}, "Mass binning"};
  ConfigurableAxis momBins{"momBins", {500, 0.f, 10.f}, "momentum binning"};
  ConfigurableAxis ResBins{"ResBins", {200, -0.2f, 0.2f}, "residual binning"};

  Configurable<bool> collSelection{"collSelection", true, "collision selection"};

  HistogramRegistry hist{"Histograms"};

  int LambdaPDG = 3122;
  int AntiLambdaPDG = -3122;
  int ProtonPDG = 2212;
  int AntiProtonPDG = -2212;
  int PosPionPDG = 211;
  int NegPionPDG = -211;

  void init(InitContext const&)
  {
    const AxisSpec mAxis{mBins, "#it{m}(p#pi) (GeV/#it{c}^{2})"};
    const AxisSpec PxAxis{momBins, "#it{p}_{x} (GeV/#it{c})"};
    const AxisSpec PyAxis{momBins, "#it{p}_{y} (GeV/#it{c})"};
    const AxisSpec PzAxis{momBins, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec PtAxis{momBins, "#it{p}_{t} (GeV/#it{c})"};
    const AxisSpec PxResidualAxis{ResBins, "(#it{p}_{x}^{rec} - #it{p}_{x}^{MC}) / #sigma_{#it{p}_{x}}"};
    const AxisSpec PyResidualAxis{ResBins, "(#it{p}_{y}^{rec} - #it{p}_{y}^{MC}) / #sigma_{#it{p}_{y}}"};
    const AxisSpec PzResidualAxis{ResBins, "(#it{p}_{z}^{rec} - #it{p}_{z}^{MC}) / #sigma_{#it{p}_{z}}"};
    const AxisSpec PtResidualAxis{ResBins, "(#it{p}_{t}^{rec} - #it{p}_{t}^{MC}) / #sigma_{#it{p}_{t}}"};

    // Lambda histograms
    hist.add("Lambda/hMassLambda", "Invariant mass Lambda", {HistType::kTH2F, {PtAxis, mAxis}});
    hist.add("Lambda/hPxResolutionProton", "#it{p}_{x} proton at V0 vtx", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("Lambda/hPxResolutionPion", "#it{p}_{x} pion at V0 vtx", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("Lambda/hPxResolutionProtonIU", "#it{p}_{x} proton at IU", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("Lambda/hPxResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("Lambda/hPyResolutionProton", "#it{p}_{y} proton at V0 vtx", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("Lambda/hPyResolutionPion", "#it{p}_{y} pion at V0 vtx", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("Lambda/hPyResolutionProtonIU", "#it{p}_{y} proton at IU", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("Lambda/hPyResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("Lambda/hPzResolutionProton", "#it{p}_{z} proton at V0 vtx", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("Lambda/hPzResolutionPion", "#it{p}_{z} pion at V0 vtx", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("Lambda/hPzResolutionProtonIU", "#it{p}_{z} proton at IU", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("Lambda/hPzResolutionPionIU", "#it{p}_{z} pion at IU", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("Lambda/hPtResolutionProton", "#it{p}_{t} proton at V0 vtx", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("Lambda/hPtResolutionPion", "#it{p}_{t} pion at V0 vtx", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("Lambda/hPtResolutionProtonIU", "#it{p}_{t} proton at IU", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("Lambda/hPtResolutionPionIU", "#it{p}_{t} pion at IU", {HistType::kTH2F, {PtAxis, PtResidualAxis}});

    // Anti-Lambda histograms
    hist.add("AntiLambda/hMassAntiLambda", "Invariant mass Anti-Lambda", {HistType::kTH2F, {PtAxis, mAxis}});
    hist.add("AntiLambda/hPxResolutionProton", "#it{p}_{x} anti-proton at V0 vtx", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("AntiLambda/hPxResolutionPion", "#it{p}_{x} pion at V0 vtx", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("AntiLambda/hPxResolutionProtonIU", "#it{p}_{x} anti-proton at IU", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("AntiLambda/hPxResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PxAxis, PxResidualAxis}});
    hist.add("AntiLambda/hPyResolutionProton", "#it{p}_{y} anti-proton at V0 vtx", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("AntiLambda/hPyResolutionPion", "#it{p}_{y} pion at V0 vtx", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("AntiLambda/hPyResolutionProtonIU", "#it{p}_{y} anti-proton at IU", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("AntiLambda/hPyResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PyAxis, PyResidualAxis}});
    hist.add("AntiLambda/hPzResolutionProton", "#it{p}_{z} anti-proton at V0 vtx", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("AntiLambda/hPzResolutionPion", "#it{p}_{z} pion at V0 vtx", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("AntiLambda/hPzResolutionProtonIU", "#it{p}_{z} anti-proton at IU", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("AntiLambda/hPzResolutionPionIU", "#it{p}_{z} pion at IU", {HistType::kTH2F, {PzAxis, PzResidualAxis}});
    hist.add("AntiLambda/hPtResolutionProton", "#it{p}_{t} anti-proton at V0 vtx", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("AntiLambda/hPtResolutionPion", "#it{p}_{t} pion at V0 vtx", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("AntiLambda/hPtResolutionProtonIU", "#it{p}_{t} anti-proton at IU", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
    hist.add("AntiLambda/hPtResolutionPionIU", "#it{p}_{t} pion at IU", {HistType::kTH2F, {PtAxis, PtResidualAxis}});
  }

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, V0DatasLabeled const& V0Datas, aod::McParticles const& mcparticles, soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels> const&)
  { 
    double PtProton, PtPion, sigmaPtProton, sigmaPtPion;
    o2::track::TrackParCov protonTrackParCov, pionTrackParCov;
    std::array<float, 21> protoncv, pioncv;

    // selection
    if (collSelection && (!collision.sel8() || abs(collision.posZ()) >= 10.))
      return;

    for (auto& v0data : V0Datas) {

      // get reconstructed info
      // double massLambda = v0data.mLambda();
      // double massAntiLambda = v0data.mAntiLambda();

      // get MC info
      if (v0data.has_mcParticle() && v0data.mcParticleId() > -1 && v0data.mcParticleId() <= mcparticles.size()) {
        auto MCv0 = v0data.mcParticle_as<aod::McParticles>();

        if (!MCv0.has_daughters()) 
          continue;

        // Lambda
        if (MCv0.pdgCode() == LambdaPDG) {
          LOG(debug) << "V0 is a Lambda.";
          const auto& protonTrackIU = v0data.posTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels>>();
          const auto& pionTrackIU = v0data.negTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels>>();

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == ProtonPDG && MCpion.pdgCode() == NegPionPDG) {

              // lambda mass
              hist.fill(HIST("Lambda/hMassLambda"), v0data.pt(), v0data.mLambda());
              // momentum resolution at Lambda vertex
              hist.fill(HIST("Lambda/hPxResolutionProton"), MCproton.px(), (v0data.pxpos() - MCproton.px()) / sqrt(v0data.covmatposdau()[9]));
              hist.fill(HIST("Lambda/hPxResolutionPion"), MCpion.px(), (v0data.pxneg() - MCpion.px()) / sqrt(v0data.covmatnegdau()[9]));
              hist.fill(HIST("Lambda/hPyResolutionProton"), MCproton.py(), (v0data.pypos() - MCproton.py()) / sqrt(v0data.covmatposdau()[14]));
              hist.fill(HIST("Lambda/hPyResolutionPion"), MCpion.py(), (v0data.pyneg() - MCpion.py()) / sqrt(v0data.covmatnegdau()[14]));
              hist.fill(HIST("Lambda/hPzResolutionProton"), MCproton.pz() ,(v0data.pzpos() - MCproton.pz()) / sqrt(v0data.covmatposdau()[20]));
              hist.fill(HIST("Lambda/hPzResolutionPion"), MCpion.pz(), (v0data.pzneg() - MCpion.pz()) / sqrt(v0data.covmatnegdau()[20]));
              // pT
              PtProton = sqrt(v0data.pxpos()*v0data.pxpos() + v0data.pypos()*v0data.pypos());
              PtPion = sqrt(v0data.pxneg()*v0data.pxneg() + v0data.pyneg()*v0data.pyneg());
              sigmaPtProton = sqrt(v0data.covmatposdau()[9] + v0data.covmatposdau()[14]);
              sigmaPtPion = sqrt(v0data.covmatnegdau()[9] + v0data.covmatnegdau()[14]);
              hist.fill(HIST("Lambda/hPtResolutionProton"), MCproton.pt(), (PtProton - MCproton.pt()) / sigmaPtProton);
              hist.fill(HIST("Lambda/hPtResolutionPion"), MCpion.pt(), (PtPion - MCpion.pt()) / sigmaPtPion);

              // momentum resolution at IU
              protonTrackParCov = getTrackParCov(protonTrackIU);
              pionTrackParCov = getTrackParCov(pionTrackIU);
              protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
              pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
              hist.fill(HIST("Lambda/hPxResolutionProtonIU"), MCproton.px(), (protonTrackIU.px() - MCproton.px()) / sqrt(protoncv[9]));
              hist.fill(HIST("Lambda/hPxResolutionPionIU"), MCpion.px(), (pionTrackIU.px() - MCpion.px()) / sqrt(pioncv[9]));
              hist.fill(HIST("Lambda/hPyResolutionProtonIU"), MCproton.py(), (protonTrackIU.py() - MCproton.py()) / sqrt(protoncv[14]));
              hist.fill(HIST("Lambda/hPyResolutionPionIU"), MCpion.py(), (pionTrackIU.py() - MCpion.py()) / sqrt(pioncv[14]));
              hist.fill(HIST("Lambda/hPzResolutionProtonIU"), MCproton.pz(), (protonTrackIU.pz() - MCproton.pz()) / sqrt(protoncv[20]));
              hist.fill(HIST("Lambda/hPzResolutionPionIU"), MCpion.pz(), (pionTrackIU.pz() - MCpion.pz()) / sqrt(pioncv[20]));
              // pT
              hist.fill(HIST("Lambda/hPtResolutionProtonIU"), MCproton.pt(), (protonTrackIU.pt() - MCproton.pt()) / sqrt(protoncv[9] + protoncv[14]));
              hist.fill(HIST("Lambda/hPtResolutionPionIU"), MCpion.pt(), (pionTrackIU.pt() - MCpion.pt()) / sqrt(pioncv[9] + pioncv[14]));
            }
          }
        } // end Lambda

        // Anti-Lambda
        if (MCv0.pdgCode() == AntiLambdaPDG) {
          LOG(debug) << "V0 is an Anti-Lambda.";
          const auto& protonTrackIU = v0data.negTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels>>();
          const auto& pionTrackIU = v0data.posTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels>>();

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == AntiProtonPDG && MCpion.pdgCode() == PosPionPDG) {

              // anti-lambda mass
              hist.fill(HIST("AntiLambda/hMassAntiLambda"), v0data.pt(), v0data.mAntiLambda());
              // momentum resolution at Lambda vertex
              hist.fill(HIST("AntiLambda/hPxResolutionProton"), MCproton.px(), (v0data.pxneg() - MCproton.px()) / sqrt(v0data.covmatnegdau()[9]));
              hist.fill(HIST("AntiLambda/hPxResolutionPion"), MCpion.px(), (v0data.pxpos() - MCpion.px()) / sqrt(v0data.covmatposdau()[9]));
              hist.fill(HIST("AntiLambda/hPyResolutionProton"), MCproton.py(), (v0data.pyneg() - MCproton.py()) / sqrt(v0data.covmatnegdau()[14]));
              hist.fill(HIST("AntiLambda/hPyResolutionPion"), MCpion.py(), (v0data.pypos() - MCpion.py()) / sqrt(v0data.covmatposdau()[14]));
              hist.fill(HIST("AntiLambda/hPzResolutionProton"), MCproton.pz(), (v0data.pzneg() - MCproton.pz()) / sqrt(v0data.covmatnegdau()[20]));
              hist.fill(HIST("AntiLambda/hPzResolutionPion"), MCpion.pz(), (v0data.pzpos() - MCpion.pz()) / sqrt(v0data.covmatposdau()[20]));
              // pT
              PtProton = sqrt(v0data.pxneg()*v0data.pxneg() + v0data.pyneg()*v0data.pyneg());
              PtPion = sqrt(v0data.pxpos()*v0data.pxpos() + v0data.pypos()*v0data.pypos());
              sigmaPtProton = sqrt(v0data.covmatnegdau()[9] + v0data.covmatnegdau()[14]);
              sigmaPtPion = sqrt(v0data.covmatposdau()[9] + v0data.covmatposdau()[14]);
              hist.fill(HIST("AntiLambda/hPtResolutionProton"), MCproton.pt(), (PtProton - MCproton.pt()) / sigmaPtProton);
              hist.fill(HIST("AntiLambda/hPtResolutionPion"), MCpion.pt(), (PtPion - MCpion.pt()) / sigmaPtPion);

              // momentum resolution at IU
              protonTrackParCov = getTrackParCov(protonTrackIU);
              pionTrackParCov = getTrackParCov(pionTrackIU);
              protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
              pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
              hist.fill(HIST("AntiLambda/hPxResolutionProtonIU"), MCproton.px(), (protonTrackIU.px() - MCproton.px()) / sqrt(protoncv[9]));
              hist.fill(HIST("AntiLambda/hPxResolutionPionIU"), MCpion.px(), (pionTrackIU.px() - MCpion.px()) / sqrt(pioncv[9]));
              hist.fill(HIST("AntiLambda/hPyResolutionProtonIU"), MCproton.py(), (protonTrackIU.py() - MCproton.py()) / sqrt(protoncv[14]));
              hist.fill(HIST("AntiLambda/hPyResolutionPionIU"), MCpion.py(), (pionTrackIU.py() - MCpion.py()) / sqrt(pioncv[14]));
              hist.fill(HIST("AntiLambda/hPzResolutionProtonIU"), MCproton.pz(), (protonTrackIU.pz() - MCproton.pz()) / sqrt(protoncv[20]));
              hist.fill(HIST("AntiLambda/hPzResolutionPionIU"), MCpion.pz(), (pionTrackIU.pz() - MCpion.pz()) / sqrt(pioncv[20]));
              // pT
              hist.fill(HIST("AntiLambda/hPtResolutionProtonIU"), MCproton.pt(), (protonTrackIU.pt() - MCproton.pt()) / sqrt(protoncv[9] + protoncv[14]));
              hist.fill(HIST("AntiLambda/hPtResolutionPionIU"), MCpion.pt(), (pionTrackIU.pt() - MCpion.pt()) / sqrt(pioncv[9] + pioncv[14]));
            }
          }
        } // end Anti-Lambda
      } // end MC
    } // end V0 loop
  }
  PROCESS_SWITCH(qaLamMomResolution, processMC, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaLamMomResolution>(cfgc)};
}
