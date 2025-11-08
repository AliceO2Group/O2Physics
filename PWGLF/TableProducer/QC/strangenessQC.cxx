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
/// \brief V0 and Cascade derived table producer for strangeness in pb-pb QC task
///
/// In case of questions please write to:
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/QC/strangenessTablesQC.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct strangenessQC {
  // Tables to produce
  Produces<aod::CascadesQC> cascadesQC;
  Produces<aod::VZerosQC> vZerosQC;

  // Histogram registries
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rGeneral{"generalInfo", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", true, "Apply sel8 event selection"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> v0setting_cospa{"v0setting_cospa", 0.97, "V0 CosPA"}; // should be double!
  Configurable<float> v0setting_radius{"v0setting_radius", 1, "v0radius"};
  Configurable<float> v0setting_rapidity{"v0setting_rapidity", 0.5, "rapidity"};

  // Configurable parameters for cascade selection
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"};        // should be double!
  Configurable<double> cascadesetting_v0cospa{"cascadesetting_v0cospa", 0.98, "Casc V0 CosPA"}; // should be double!
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "DCA cascade daughters"};
  Configurable<float> cascadesetting_dcav0dau{"cascadesetting_dcav0dau", 1.0, "DCA Cascade's V0 Daughters"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "DCA bachelor to PV"};
  Configurable<float> cascadesetting_dcapostopv{"cascadesetting_dcapostopv", 0.2, "DCA Casc. V0's pos to PV"};
  Configurable<float> cascadesetting_dcanegtopv{"cascadesetting_dcanegtopv", 0.2, "DCA Casc V0's neg to PV"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "minimum V0 DCA to PV"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 1.0, "cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "v0 mass window"};
  Configurable<float> cascadesetting_v0radius{"cascadesetting_v0radius", 0.9, "v0 radius"};
  Configurable<float> cascadesetting_rapidity{"cascadesetting_rapidity", 0.5, "rapidity"};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {100, -15., 15., "vrtx_{Z} [cm]"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // Cut summary
    rGeneral.add("selectionSummary", "selectionSummary", HistType::kTH1F, {{19, -0.5, 18.5}});
    TString CutLabelSummary[19] = {"cutzvertex", "v0_dcav0dau", "v0_dcapostopv", "v0_dcanegtopv", "v0_cospa", "v0_radius", "v0_rapidity",
                                   "casc_cospa", "casc_v0cospa", "casc_dcacascdau", "casc_dcav0dau", "casc_dcabachtopv", "casc_dcapostopv", "casc_dcanegtopv", "casc_mindcav0topv", "casc_cascradius", "casc_v0masswindow", "casc_v0radius", "casc_rapidity"};
    for (Int_t n = 1; n <= rGeneral.get<TH1>(HIST("selectionSummary"))->GetNbinsX(); n++) {
      rGeneral.get<TH1>(HIST("selectionSummary"))->GetXaxis()->SetBinLabel(n, CutLabelSummary[n - 1]);
    }
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(1, cutzvertex);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(2, v0setting_dcav0dau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(3, v0setting_dcapostopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(4, v0setting_dcanegtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(5, v0setting_cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(6, v0setting_radius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(7, v0setting_rapidity);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(8, cascadesetting_cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(9, cascadesetting_v0cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(10, cascadesetting_dcacascdau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(11, cascadesetting_dcav0dau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(12, cascadesetting_dcabachtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(13, cascadesetting_dcapostopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(14, cascadesetting_dcanegtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(15, cascadesetting_mindcav0topv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(16, cascadesetting_cascradius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(17, cascadesetting_v0masswindow);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(18, cascadesetting_v0radius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(19, cascadesetting_rapidity);
  }

  // Event selection
  Filter eventFilter = (sel8 && o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau &&
                        aod::v0data::v0cosPA > v0setting_cospa);

  // Filters on Cascades
  Filter preFilterCascades = (nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascdata::dcaV0daughters < cascadesetting_dcav0dau &&
                              nabs(aod::cascdata::dcapostopv) > cascadesetting_dcapostopv &&
                              nabs(aod::cascdata::dcanegtopv) > cascadesetting_dcanegtopv &&
                              aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   soa::Filtered<aod::CascDatas> const& Cascades,
                   soa::Filtered<aod::V0Datas> const& V0s,
                   aod::V0Datas const&, // it's needed to access the full table of V0s (not the filtered one) to make sure all the V0s related to cascades are present
                   aod::V0sLinked const&,
                   DaughterTracks const&)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // v0-s loop
    for (const auto& v0 : V0s) {
      // Decay radius check
      if (v0.v0radius() < v0setting_radius) {
        continue;
      }
      // Rapidity check
      if (TMath::Abs(v0.yK0Short()) > v0setting_rapidity &&
          TMath::Abs(v0.yLambda()) > v0setting_rapidity) {
        continue;
      }

      const auto& posDau = v0.posTrack_as<DaughterTracks>();
      const auto& negDau = v0.negTrack_as<DaughterTracks>();

      int posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (posDau.itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (negDau.itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      float decayLength = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(v0.px(), v0.py(), v0.pz());
      float cTauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * pdgDB->Mass(3122);
      float cTauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * pdgDB->Mass(310);

      vZerosQC(v0.yK0Short(), v0.yLambda(),
               v0.dcaV0daughters(), v0.dcav0topv(),
               v0.dcapostopv(), v0.dcanegtopv(),
               v0.v0radius(), v0.v0cosPA(),
               posDau.tpcNSigmaPi(), posDau.tpcNSigmaPr(),
               negDau.tpcNSigmaPi(), negDau.tpcNSigmaPr(),
               decayLength, cTauLambda, cTauK0s,
               posITSNhits, negITSNhits,
               v0.pt(), v0.eta(),
               v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(),
               posDau.eta(), negDau.eta(),
               posDau.phi(), negDau.phi());
    }

    // cascades loop
    for (const auto& casc : Cascades) {
      // Rapidity check
      if (TMath::Abs(casc.yXi()) > cascadesetting_rapidity &&
          TMath::Abs(casc.yOmega()) > cascadesetting_rapidity) {
        continue;
      }

      float casc_v0cospa = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float casc_cospa = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float casc_dcav0topv = casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ());

      // Cut on dynamic columns
      if (casc.v0radius() < cascadesetting_v0radius ||
          casc.cascradius() < cascadesetting_cascradius ||
          casc_v0cospa < cascadesetting_v0cospa ||
          casc_cospa < cascadesetting_cospa ||
          TMath::Abs(casc_dcav0topv) < cascadesetting_mindcav0topv) {
        continue;
      }

      // Cut on v0 invariant mass
      if (TMath::Abs(casc.mLambda() - pdgDB->Mass(3122)) > cascadesetting_v0masswindow)
        continue;

      const auto& bachDau = casc.bachelor_as<DaughterTracks>();
      const auto& posDau = casc.posTrack_as<DaughterTracks>();
      const auto& negDau = casc.negTrack_as<DaughterTracks>();

      float cascDecayLength = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascTotalMomentum = std::hypot(casc.px(), casc.py(), casc.pz());
      float CtauXi = cascDecayLength / (cascTotalMomentum + 1e-13) * pdgDB->Mass(3312);
      float CtauOmega = cascDecayLength / (cascTotalMomentum + 1e-13) * pdgDB->Mass(3334);

      float v0DecayLength = std::hypot(casc.xlambda() - casc.x(), casc.ylambda() - casc.y(), casc.zlambda() - casc.z());
      float v0TotalMomentum = std::hypot(casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg());
      float CtauV0 = v0DecayLength / (v0TotalMomentum + 1e-13) * pdgDB->Mass(3122);

      cascadesQC(casc.sign(), casc.yXi(), casc.yOmega(),
                 casc_cospa, casc_v0cospa,
                 casc.cascradius(), casc.v0radius(),
                 cascDecayLength, CtauXi, CtauOmega, CtauV0,
                 casc.dcaV0daughters(), casc.dcacascdaughters(), casc_dcav0topv,
                 casc.dcabachtopv(), casc.dcapostopv(), casc.dcanegtopv(),
                 posDau.tpcNSigmaPi(), posDau.tpcNSigmaPr(),
                 negDau.tpcNSigmaPi(), negDau.tpcNSigmaPr(),
                 bachDau.tpcNSigmaPi(), bachDau.tpcNSigmaKa(),
                 casc.pt(), casc.eta(),
                 casc.mLambda(), casc.mOmega(), casc.mXi());
    }
  }
  PROCESS_SWITCH(strangenessQC, processData, "Process Run 3 data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessQC>(cfgc, TaskName{"lf-strangenessqc"})};
}
