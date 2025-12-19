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
// This task is designed to do QA to the TOF PID applied to strangeness
// in the regular framework

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

struct strangepidqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "V0 radius (cm)"};

  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 50.0f, 70.0f, 100.0f}, "FT0C centrality"};

  // Invariant Mass
  ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.497f - 0.050f, 0.497f + 0.050f}, "M_{K0s} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.08f, 1.16f}, "M_{#Lambda} (GeV/c^{2})"};
  AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

  // Length axis
  ConfigurableAxis axisLength{"axisLength", {600, 0.0f, +600.0f}, "track Length (cm)"};

  // TOF cut axis
  ConfigurableAxis axisTOFCut{"axisTOFCut", {100, 0.0f, +10000.0f}, "TOF compat. cut (ps)"};

  // nsigma axis
  ConfigurableAxis axisNSigma{"axisNSigma", {60, -3.0f, +3.0f}, "NSigma"};

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Selection criteria for cascade analysis
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.95, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Track configurables
  Configurable<int> tpcCrossedRows{"tpcCrossedRows", 70, "minimum number of TPC rows requirement"};
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 4, "TPC NSigma bachelor (>10 is no cut)"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 4, "TPC NSigma proton <- lambda (>10 is no cut)"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 4, "TPC NSigma pion <- lambda (>10 is no cut)"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Direct test configurables
  Configurable<float> massWindowForNSigmaPlots{"massWindowForNSigmaPlots", 0.005f, "mass window for Nsigma comparison plots"};
  Configurable<float> tofNsigmaCompatibility{"tofNsigmaCompatibility", 4, "compatibility check for V0s"};
  Configurable<float> tofNsigmaCompatibilityCascades{"tofNsigmaCompatibilityCascades", 4, "compatibility check for cascades"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, 100.0f}});

    histos.add("h1dMassK0Short", "h1dMassK0Short", {HistType::kTH1F, {axisK0ShortMass}});
    histos.add("h1dMassLambda", "h1dMassLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassAntiLambda", "h1dMassAntiLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassCompatibleK0Short", "h1dMassCompatibleK0Short", {HistType::kTH1F, {axisK0ShortMass}});
    histos.add("h1dMassCompatibleLambda", "h1dMassCompatibleLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassCompatibleAntiLambda", "h1dMassCompatibleAntiLambda", {HistType::kTH1F, {axisLambdaMass}});

    histos.add("h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {centAxis, axisPt, axisK0ShortMass}});
    histos.add("h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dMassCompatibleK0Short", "h3dMassCompatibleK0Short", {HistType::kTH3F, {centAxis, axisPt, axisK0ShortMass}});
    histos.add("h3dMassCompatibleLambda", "h3dMassCompatibleLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dMassCompatibleAntiLambda", "h3dMassCompatibleAntiLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});

    // cross-check if compatibility requested with primary TOF instead (requires non-derived)
    histos.add("h3dPrimaryTOFMassCompatibleK0Short", "h3dPrimaryTOFMassCompatibleK0Short", {HistType::kTH3F, {centAxis, axisPt, axisK0ShortMass}});
    histos.add("h3dPrimaryTOFMassCompatibleLambda", "h3dPrimaryTOFMassCompatibleLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dPrimaryTOFMassCompatibleAntiLambda", "h3dPrimaryTOFMassCompatibleAntiLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});

    // plot Nsigma: primary vs secondary TOF vs pT (use narrow window around mass)
    histos.add("h3dNSigmasLaPr", "h3dNSigmasLaPr", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
    histos.add("h3dNSigmasLaPi", "h3dNSigmasLaPi", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
    histos.add("h3dNSigmasK0Pi", "h3dNSigmasK0Pi", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});

    if (doprocessCascades || doprocessCascadesNonDerived) {
      histos.add("h1dMassXiMinus", "h1dMassXiMinus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassXiPlus", "h1dMassXiPlus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassOmegaMinus", "h1dMassOmegaMinus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassOmegaPlus", "h1dMassOmegaPlus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassCompatibleXiMinus", "h1dMassCompatibleXiMinus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassCompatibleXiPlus", "h1dMassCompatibleXiPlus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassCompatibleOmegaMinus", "h1dMassCompatibleOmegaMinus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassCompatibleOmegaPlus", "h1dMassCompatibleOmegaPlus", {HistType::kTH1F, {massAxisOmega}});

      histos.add("h3dMassXiMinus", "h3dMassXiMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassXiPlus", "h3dMassXiPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dMassCompatibleXiMinus", "h3dMassCompatibleXiMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassCompatibleXiPlus", "h3dMassCompatibleXiPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassCompatibleOmegaMinus", "h3dMassCompatibleOmegaMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dMassCompatibleOmegaPlus", "h3dMassCompatibleOmegaPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});

      // cross-check if compatibility requested with primary TOF instead (requires non-derived)
      histos.add("h3dPrimaryTOFMassCompatibleXiMinus", "h3dPrimaryTOFMassCompatibleXiMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dPrimaryTOFMassCompatibleXiPlus", "h3dPrimaryTOFMassCompatibleXiPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dPrimaryTOFMassCompatibleOmegaMinus", "h3dPrimaryTOFMassCompatibleOmegaMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dPrimaryTOFMassCompatibleOmegaPlus", "h3dPrimaryTOFMassCompatibleOmegaPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});

      // plot Nsigma: primary vs secondary TOF vs pT (use narrow window around mass)
      histos.add("h3dNSigmasXiLaPr", "h3dNSigmasXiLaPr", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
      histos.add("h3dNSigmasXiLaPi", "h3dNSigmasXiLaPi", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
      histos.add("h3dNSigmasXiPi", "h3dNSigmasXiPi", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});

      histos.add("h3dNSigmasOmLaPr", "h3dNSigmasOmLaPr", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
      histos.add("h3dNSigmasOmLaPi", "h3dNSigmasOmLaPi", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
      histos.add("h3dNSigmasOmKa", "h3dNSigmasOmKa", {HistType::kTH3F, {axisNSigma, axisNSigma, axisPt}});
    }
  }

  void processReal(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFBetas, aod::V0TOFDebugs, aod::V0TOFNSigmas> const& v0s, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    for (auto& lambda : v0s) {

      if (TMath::Abs(lambda.eta()) > 0.5)
        continue;

      auto negExtra = lambda.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = lambda.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
        histos.fill(HIST("h1dMassLambda"), lambda.mLambda());
        if (lambda.tofLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
          histos.fill(HIST("h1dMassCompatibleLambda"), lambda.mLambda());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
        histos.fill(HIST("h1dMassAntiLambda"), lambda.mAntiLambda());
        if (lambda.tofAntiLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
          histos.fill(HIST("h1dMassCompatibleAntiLambda"), lambda.mAntiLambda());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
        histos.fill(HIST("h1dMassK0Short"), lambda.mK0Short());
        if (lambda.tofK0ShortCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
          histos.fill(HIST("h1dMassCompatibleK0Short"), lambda.mK0Short());
        }
      }
    }
  }

  // ____________________________________________________________________________
  // QA TOF NSigma quantities
  Filter preFilter =
    nabs(aod::cascdata::dcapostopv) > v0setting_dcapostopv&& nabs(aod::cascdata::dcanegtopv) > v0setting_dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < v0setting_dcav0dau&& aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau;

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents> const& collisions, soa::Filtered<soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascTOFNSigmas>> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    for (auto& casc : Cascades) {
      auto col = collisions.rawIteratorAt(casc.straCollisionId());

      // major selections here
      if (casc.v0radius() > v0setting_radius &&
          casc.cascradius() > cascadesetting_cascradius &&
          casc.v0cosPA(col.posX(), col.posY(), col.posZ()) > v0setting_cospa &&
          casc.casccosPA(col.posX(), col.posY(), col.posZ()) > cascadesetting_cospa &&
          casc.dcav0topv(col.posX(), col.posY(), col.posZ()) > cascadesetting_mindcav0topv &&
          TMath::Abs(casc.mLambda() - 1.115683) < cascadesetting_v0masswindow) {

        auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
        auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
        auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

        if (negExtra.tpcCrossedRows() < tpcCrossedRows || posExtra.tpcCrossedRows() < tpcCrossedRows || bachExtra.tpcCrossedRows() < tpcCrossedRows)
          continue;

        if (casc.sign() < 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiMinus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiMinus"), casc.mXi());
            }
          }
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaMinus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaMinus"), casc.mOmega());
            }
          }
        } else {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiPlus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiPlus"), casc.mXi());
            }
          }

          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaPlus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaPlus"), casc.mOmega());
            }
          }
        }
      }
    }
  }

  void processRealNonDerived(soa::Join<aod::Collisions, aod::CentFT0Cs> const& collisions, soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFBetas, aod::V0TOFDebugs, aod::V0TOFNSigmas> const& v0s, soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr> const&)
  {
    for (auto const& col : collisions) {
      histos.fill(HIST("hEventCentrality"), col.centFT0C());
    }

    for (auto& lambda : v0s) {
      auto coll = collisions.rawIteratorAt(lambda.collisionId());

      if (TMath::Abs(lambda.eta()) > 0.5)
        continue;

      auto negExtra = lambda.negTrack_as<soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>>();
      auto posExtra = lambda.posTrack_as<soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>>();

      bool primaryTOFcompatible_Lambda =
        (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value)) &&
        (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value));

      bool primaryTOFcompatible_AntiLambda =
        (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
        (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value));

      bool primaryTOFcompatible_K0Short =
        (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
        (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value));

      if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
        histos.fill(HIST("h1dMassLambda"), lambda.mLambda());
        if (lambda.tofLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
          histos.fill(HIST("h1dMassCompatibleLambda"), lambda.mLambda());
        }
        if (primaryTOFcompatible_Lambda) {
          histos.fill(HIST("h3dPrimaryTOFMassCompatibleLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
        }
        if (std::abs(lambda.mLambda() - o2::constants::physics::MassLambda) < massWindowForNSigmaPlots.value) {
          histos.fill(HIST("h3dNSigmasLaPr"), lambda.tofNSigmaLaPr(), posExtra.tofNSigmaPr(), lambda.pt());
          histos.fill(HIST("h3dNSigmasLaPi"), lambda.tofNSigmaLaPi(), negExtra.tofNSigmaPi(), lambda.pt());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
        histos.fill(HIST("h1dMassAntiLambda"), lambda.mAntiLambda());
        if (lambda.tofAntiLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
          histos.fill(HIST("h1dMassCompatibleAntiLambda"), lambda.mAntiLambda());
        }
        if (primaryTOFcompatible_AntiLambda) {
          histos.fill(HIST("h3dPrimaryTOFMassCompatibleAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
        histos.fill(HIST("h1dMassK0Short"), lambda.mK0Short());
        if (lambda.tofK0ShortCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
          histos.fill(HIST("h1dMassCompatibleK0Short"), lambda.mK0Short());
        }
        if (primaryTOFcompatible_K0Short) {
          histos.fill(HIST("h3dPrimaryTOFMassCompatibleK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
        }
        if (std::abs(lambda.mK0Short() - o2::constants::physics::MassK0Short) < massWindowForNSigmaPlots.value) {
          histos.fill(HIST("h3dNSigmasK0Pi"), lambda.tofNSigmaLaPi(), posExtra.tofNSigmaPi(), lambda.pt());
          histos.fill(HIST("h3dNSigmasK0Pi"), lambda.tofNSigmaLaPi(), negExtra.tofNSigmaPi(), lambda.pt());
        }
      }
    }
  }

  // to test original data (not derived) as well, compare with primary TOF Nsigmas
  // don't do grouping, faster to simply stream through
  void processCascadesNonDerived(soa::Join<aod::Collisions, aod::CentFT0Cs> const& collisions, soa::Filtered<soa::Join<aod::CascDatas, aod::CascTOFNSigmas>> const& Cascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr> const&)
  {
    for (auto& casc : Cascades) {
      auto col = collisions.rawIteratorAt(casc.collisionId());

      // major selections here
      if (casc.v0radius() > v0setting_radius &&
          casc.cascradius() > cascadesetting_cascradius &&
          casc.v0cosPA(col.posX(), col.posY(), col.posZ()) > v0setting_cospa &&
          casc.casccosPA(col.posX(), col.posY(), col.posZ()) > cascadesetting_cospa &&
          casc.dcav0topv(col.posX(), col.posY(), col.posZ()) > cascadesetting_mindcav0topv &&
          TMath::Abs(casc.mLambda() - 1.115683) < cascadesetting_v0masswindow) {

        auto negExtra = casc.negTrack_as<soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>>();
        auto posExtra = casc.posTrack_as<soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>>();
        auto bachExtra = casc.bachelor_as<soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>>();

        if (negExtra.tpcNClsCrossedRows() < tpcCrossedRows || posExtra.tpcNClsCrossedRows() < tpcCrossedRows || bachExtra.tpcNClsCrossedRows() < tpcCrossedRows) {
          continue;
        }

        bool primaryTOFcompatible_XiMinus =
          (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value)) &&
          (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
          (!bachExtra.hasTOF() || (std::fabs(bachExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value));

        bool primaryTOFcompatible_XiPlus =
          (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
          (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value)) &&
          (!bachExtra.hasTOF() || (std::fabs(bachExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value));

        bool primaryTOFcompatible_OmegaMinus =
          (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value)) &&
          (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
          (!bachExtra.hasTOF() || (std::fabs(bachExtra.tofNSigmaKa()) < tofNsigmaCompatibilityCascades.value));

        bool primaryTOFcompatible_OmegaPlus =
          (!posExtra.hasTOF() || (std::fabs(posExtra.tofNSigmaPi()) < tofNsigmaCompatibilityCascades.value)) &&
          (!negExtra.hasTOF() || (std::fabs(negExtra.tofNSigmaPr()) < tofNsigmaCompatibilityCascades.value)) &&
          (!bachExtra.hasTOF() || (std::fabs(bachExtra.tofNSigmaKa()) < tofNsigmaCompatibilityCascades.value));

        if (casc.sign() < 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiMinus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiMinus"), casc.mXi());
            }
            if (primaryTOFcompatible_XiMinus) {
              histos.fill(HIST("h3dPrimaryTOFMassCompatibleXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
            }
            if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) < massWindowForNSigmaPlots.value) {
              histos.fill(HIST("h3dNSigmasXiLaPr"), casc.tofNSigmaXiLaPr(), posExtra.tofNSigmaPr(), casc.pt());
              histos.fill(HIST("h3dNSigmasXiLaPi"), casc.tofNSigmaXiLaPi(), negExtra.tofNSigmaPi(), casc.pt());
              histos.fill(HIST("h3dNSigmasXiPi"), casc.tofNSigmaXiPi(), bachExtra.tofNSigmaPi(), casc.pt());
            }
          }
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaMinus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaMinus"), casc.mOmega());
            }
            if (primaryTOFcompatible_OmegaMinus) {
              histos.fill(HIST("h3dPrimaryTOFMassCompatibleOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
            }
            if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < massWindowForNSigmaPlots.value) {
              histos.fill(HIST("h3dNSigmasOmLaPr"), casc.tofNSigmaOmLaPr(), posExtra.tofNSigmaPr(), casc.pt());
              histos.fill(HIST("h3dNSigmasOmLaPi"), casc.tofNSigmaOmLaPi(), negExtra.tofNSigmaPi(), casc.pt());
              histos.fill(HIST("h3dNSigmasOmKa"), casc.tofNSigmaOmKa(), bachExtra.tofNSigmaKa(), casc.pt());
            }
          }
        } else {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiPlus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiPlus"), casc.mXi());
            }
            if (primaryTOFcompatible_XiPlus) {
              histos.fill(HIST("h3dPrimaryTOFMassCompatibleXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
            }
          }

          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaPlus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaPlus"), casc.mOmega());
            }
            if (primaryTOFcompatible_OmegaPlus) {
              histos.fill(HIST("h3dPrimaryTOFMassCompatibleOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(strangepidqa, processReal, "Produce real information", true);
  PROCESS_SWITCH(strangepidqa, processCascades, "Process real cascades", true);

  // non-derived options
  PROCESS_SWITCH(strangepidqa, processRealNonDerived, "Process real cascades from non-derived data", true);
  PROCESS_SWITCH(strangepidqa, processCascadesNonDerived, "Process real cascades from non-derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangepidqa>(cfgc)};
}
