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

  HistogramRegistry hist{"Histograms"};

  // collision selection
  Configurable<bool> collSelection{"collSelection", true, "collision selection"};
  Filter collisionFilter = (collSelection && aod::evsel::sel8 == true && nabs(aod::collision::posZ) < 10.);

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
    hist.add("Lambda/hPxResolutionProton", "#it{p}_{x} proton at V0 vtx", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("Lambda/hPxResolutionPion", "#it{p}_{x} pion at V0 vtx", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("Lambda/hPxResolutionProtonIU", "#it{p}_{x} proton at IU", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("Lambda/hPxResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("Lambda/hPyResolutionProton", "#it{p}_{y} proton at V0 vtx", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("Lambda/hPyResolutionPion", "#it{p}_{y} pion at V0 vtx", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("Lambda/hPyResolutionProtonIU", "#it{p}_{y} proton at IU", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("Lambda/hPyResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("Lambda/hPzResolutionProton", "#it{p}_{z} proton at V0 vtx", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("Lambda/hPzResolutionPion", "#it{p}_{z} pion at V0 vtx", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("Lambda/hPzResolutionProtonIU", "#it{p}_{z} proton at IU", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("Lambda/hPzResolutionPionIU", "#it{p}_{z} pion at IU", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("Lambda/hPtResolutionProton", "#it{p}_{t} proton at V0 vtx", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("Lambda/hPtResolutionPion", "#it{p}_{t} pion at V0 vtx", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("Lambda/hPtResolutionProtonIU", "#it{p}_{t} proton at IU", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("Lambda/hPtResolutionPionIU", "#it{p}_{t} pion at IU", {HistType::kTH2F, {PtResidualAxis, PtAxis}});

    // Anti-Lambda histograms
    hist.add("AntiLambda/hPxResolutionProton", "#it{p}_{x} anti-proton at V0 vtx", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("AntiLambda/hPxResolutionPion", "#it{p}_{x} pion at V0 vtx", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("AntiLambda/hPxResolutionProtonIU", "#it{p}_{x} anti-proton at IU", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("AntiLambda/hPxResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PxResidualAxis, PxAxis}});
    hist.add("AntiLambda/hPyResolutionProton", "#it{p}_{y} anti-proton at V0 vtx", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("AntiLambda/hPyResolutionPion", "#it{p}_{y} pion at V0 vtx", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("AntiLambda/hPyResolutionProtonIU", "#it{p}_{y} anti-proton at IU", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("AntiLambda/hPyResolutionPionIU", "#it{p}_{x} pion at IU", {HistType::kTH2F, {PyResidualAxis, PyAxis}});
    hist.add("AntiLambda/hPzResolutionProton", "#it{p}_{z} anti-proton at V0 vtx", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("AntiLambda/hPzResolutionPion", "#it{p}_{z} pion at V0 vtx", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("AntiLambda/hPzResolutionProtonIU", "#it{p}_{z} anti-proton at IU", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("AntiLambda/hPzResolutionPionIU", "#it{p}_{z} pion at IU", {HistType::kTH2F, {PzResidualAxis, PzAxis}});
    hist.add("AntiLambda/hPtResolutionProton", "#it{p}_{t} anti-proton at V0 vtx", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("AntiLambda/hPtResolutionPion", "#it{p}_{t} pion at V0 vtx", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("AntiLambda/hPtResolutionProtonIU", "#it{p}_{t} anti-proton at IU", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
    hist.add("AntiLambda/hPtResolutionPionIU", "#it{p}_{t} pion at IU", {HistType::kTH2F, {PtResidualAxis, PtAxis}});
  }

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, V0DatasLabeled const& V0Datas, aod::McParticles const& mcparticles, soa::Join<aod::TracksIU, aod::TracksCovIU> const&)
  { 
    double PtProton, PtPion, sigmaPtProton, sigmaPtPion;
    o2::track::TrackParCov protonTrackParCov, pionTrackParCov;
    std::array<float, 21> protoncv, pioncv;

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
          LOG(info) << "V0 is a Lambda.";
          for (auto& Daughter0 : MCv0.daughters_as<aod::McParticles>()) {
            for (auto& Daughter1 : MCv0.daughters_as<aod::McParticles>()) {
              if (Daughter0.pdgCode() == ProtonPDG && Daughter1.pdgCode() == NegPionPDG) {
                hist.fill(HIST("hMassLambda"), v0data.mLambda(), v0data.pt());
                // momentum resolution at Lambda vertex
                hist.fill(HIST("Lambda/hPxResolutionProton"), (v0data.pxpos() - Daughter0.px()) / sqrt(v0data.covmatposdau()[9]), Daughter0.px());
                hist.fill(HIST("Lambda/hPxResolutionPion"), (v0data.pxneg() - Daughter1.px()) / sqrt(v0data.covmatnegdau()[9]), Daughter1.px());
                hist.fill(HIST("Lambda/hPyResolutionProton"), (v0data.pypos() - Daughter0.py()) / sqrt(v0data.covmatposdau()[14]), Daughter0.py());
                hist.fill(HIST("Lambda/hPyResolutionPion"), (v0data.pyneg() - Daughter1.py()) / sqrt(v0data.covmatnegdau()[14]), Daughter1.py());
                hist.fill(HIST("Lambda/hPzResolutionProton"), (v0data.pzpos() - Daughter0.pz()) / sqrt(v0data.covmatposdau()[20]), Daughter0.pz());
                hist.fill(HIST("Lambda/hPzResolutionPion"), (v0data.pzneg() - Daughter1.pz()) / sqrt(v0data.covmatnegdau()[20]), Daughter1.pz());
                // pT
                PtProton = sqrt(v0data.pxpos()*v0data.pxpos() + v0data.pypos()*v0data.pypos());
                PtPion = sqrt(v0data.pxneg()*v0data.pxneg() + v0data.pyneg()*v0data.pyneg());
                sigmaPtProton = sqrt(v0data.covmatposdau()[9] + v0data.covmatposdau()[14]);
                sigmaPtPion = sqrt(v0data.covmatnegdau()[9] + v0data.covmatnegdau()[14]);
                hist.fill(HIST("Lambda/hPtResolutionProton"), (PtProton - Daughter0.pt()) / sigmaPtProton, Daughter0.pt());
                hist.fill(HIST("Lambda/hPtResolutionPion"), (PtPion - Daughter1.pt()) / sigmaPtPion, Daughter1.pt());

                // momentum resolution at IU
                auto protonTrackIU = v0data.posTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU>>();
                auto pionTrackIU = v0data.negTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU>>();
                protonTrackParCov = getTrackParCov(protonTrackIU);
                pionTrackParCov = getTrackParCov(pionTrackIU);
                protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
                pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
                hist.fill(HIST("Lambda/hPxResolutionProtonIU"), (protonTrackIU.px() - Daughter0.px()) / sqrt(protoncv[9]), Daughter0.px());
                hist.fill(HIST("Lambda/hPxResolutionPionIU"), (pionTrackIU.px() - Daughter1.px()) / sqrt(pioncv[9]), Daughter1.px());
                hist.fill(HIST("Lambda/hPyResolutionProtonIU"), (protonTrackIU.py() - Daughter0.py()) / sqrt(protoncv[14]), Daughter0.py());
                hist.fill(HIST("Lambda/hPyResolutionPionIU"), (pionTrackIU.py() - Daughter1.py()) / sqrt(pioncv[14]), Daughter1.py());
                hist.fill(HIST("Lambda/hPzResolutionProtonIU"), (protonTrackIU.pz() - Daughter0.pz()) / sqrt(protoncv[20]), Daughter0.pz());
                hist.fill(HIST("Lambda/hPzResolutionPionIU"), (pionTrackIU.pz() - Daughter1.pz()) / sqrt(pioncv[20]), Daughter1.pz());
                // pT
                hist.fill(HIST("Lambda/hPtResolutionProtonIU"), (protonTrackIU.pt() - Daughter0.pt()) / sqrt(protoncv[9] + protoncv[14]), Daughter0.pt());
                hist.fill(HIST("Lambda/hPtResolutionPionIU"), (pionTrackIU.pt() - Daughter1.pt()) / sqrt(pioncv[9] + pioncv[14]), Daughter1.pt());
              }
            } // end daughter 1
          } // end daughter 0
        } // end Lambda

        // Anti-Lambda
        if (MCv0.pdgCode() == AntiLambdaPDG) {
          LOG(info) << "V0 is an Anti-Lambda.";
          for (auto& Daughter0 : MCv0.daughters_as<aod::McParticles>()) {
            for (auto& Daughter1 : MCv0.daughters_as<aod::McParticles>()) {
              if (Daughter0.pdgCode() == AntiProtonPDG && Daughter1.pdgCode() == PosPionPDG) {
                hist.fill(HIST("hMassAntiLambda"), v0data.mAntiLambda(), v0data.pt());
                // momentum resolution at Lambda vertex
                hist.fill(HIST("AntiLambda/hPxResolutionProton"), (v0data.pxneg() - Daughter0.px()) / sqrt(v0data.covmatnegdau()[9]), Daughter0.px());
                hist.fill(HIST("AntiLambda/hPxResolutionPion"), (v0data.pxpos() - Daughter1.px()) / sqrt(v0data.covmatposdau()[9]), Daughter1.px());
                hist.fill(HIST("AntiLambda/hPyResolutionProton"), (v0data.pyneg() - Daughter0.py()) / sqrt(v0data.covmatnegdau()[14]), Daughter0.py());
                hist.fill(HIST("AntiLambda/hPyResolutionPion"), (v0data.pypos() - Daughter1.py()) / sqrt(v0data.covmatposdau()[14]), Daughter1.py());
                hist.fill(HIST("AntiLambda/hPzResolutionProton"), (v0data.pzneg() - Daughter0.pz()) / sqrt(v0data.covmatnegdau()[20]), Daughter0.pz());
                hist.fill(HIST("AntiLambda/hPzResolutionPion"), (v0data.pzpos() - Daughter1.pz()) / sqrt(v0data.covmatposdau()[20]), Daughter1.pz());
                // pT
                PtProton = sqrt(v0data.pxneg()*v0data.pxneg() + v0data.pyneg()*v0data.pyneg());
                PtPion = sqrt(v0data.pxpos()*v0data.pxpos() + v0data.pypos()*v0data.pypos());
                sigmaPtProton = sqrt(v0data.covmatnegdau()[9] + v0data.covmatnegdau()[14]);
                sigmaPtPion = sqrt(v0data.covmatposdau()[9] + v0data.covmatposdau()[14]);
                hist.fill(HIST("AntiLambda/hPtResolutionProton"), (PtProton - Daughter0.pt()) / sigmaPtProton, Daughter0.pt());
                hist.fill(HIST("AntiLambda/hPtResolutionPion"), (PtPion - Daughter1.pt()) / sigmaPtPion, Daughter1.pt());

                // momentum resolution at IU
                auto protonTrackIU = v0data.negTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU>>();
                auto pionTrackIU = v0data.posTrack_as<soa::Join<aod::TracksIU, aod::TracksCovIU>>();
                protonTrackParCov = getTrackParCov(protonTrackIU);
                pionTrackParCov = getTrackParCov(pionTrackIU);
                protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
                pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
                hist.fill(HIST("AntiLambda/hPxResolutionProtonIU"), (protonTrackIU.px() - Daughter0.px()) / sqrt(protoncv[9]), Daughter0.px());
                hist.fill(HIST("AntiLambda/hPxResolutionPionIU"), (pionTrackIU.px() - Daughter1.px()) / sqrt(pioncv[9]), Daughter1.px());
                hist.fill(HIST("AntiLambda/hPyResolutionProtonIU"), (protonTrackIU.py() - Daughter0.py()) / sqrt(protoncv[14]), Daughter0.py());
                hist.fill(HIST("AntiLambda/hPyResolutionPionIU"), (pionTrackIU.py() - Daughter1.py()) / sqrt(pioncv[14]), Daughter1.py());
                hist.fill(HIST("AntiLambda/hPzResolutionProtonIU"), (protonTrackIU.pz() - Daughter0.pz()) / sqrt(protoncv[20]), Daughter0.pz());
                hist.fill(HIST("AntiLambda/hPzResolutionPionIU"), (pionTrackIU.pz() - Daughter1.pz()) / sqrt(pioncv[20]), Daughter1.pz());
                // pT
                hist.fill(HIST("AntiLambda/hPtResolutionProtonIU"), (protonTrackIU.pt() - Daughter0.pt()) / sqrt(protoncv[9] + protoncv[14]), Daughter0.pt());
                hist.fill(HIST("AntiLambda/hPtResolutionPionIU"), (pionTrackIU.pt() - Daughter1.pt()) / sqrt(pioncv[9] + pioncv[14]), Daughter1.pt());
              }
            } // end daughter 1
          } // end daughter 0
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
