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
/// \ post processing for Cascade analysis runing on derived data
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

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

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct derivedCascadeAnalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.222f, 1.422f}, ""};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.572f, 1.772f}, ""};

  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};

  Configurable<float> minPt{"minPt", 0.0f, "minPt"};
  Configurable<float> masswin{"masswin", 0.1, "Mass window limit"};
  Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity acceptance"};
  Configurable<float> etaDauCut{"etaDauCut", 0.8, "Pseudorapidity acceptance of the cascade daughters"};
  Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.1, "DCA of baryon doughter track To PV"};
  Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.2, "DCA of meson doughter track To PV"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.1, "DCA Bach To PV"};
  Configurable<float> casCosPaPtParameter{"casCosPaPtParameter", 0.341715, "Parameter for pt dependent cos PA cut"};
  Configurable<double> casccospa{"casccospa", 0.95, "Casc CosPA"};
  Configurable<float> v0CosPaPtParameter{"v0CosPaPtParameter", 0.341715, "Parameter for pt dependent cos PA cut of the V0 daughter"};
  Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"};
  Configurable<float> dcacascdau{"dcacascdau", 1, "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> dcaV0ToPV{"dcaV0ToPV", 0.06, "DCA V0 To PV"};
  Configurable<float> minRadius{"minRadius", 1.2f, "minRadius"};
  Configurable<float> maxRadius{"maxRadius", 100.0f, "maxRadius"};
  Configurable<float> minV0Radius{"minV0Radius", 3.0f, "V0 transverse decay radius, minimum"};
  Configurable<float> maxV0Radius{"maxV0Radius", 100.0f, "V0 transverse decay radius, maximum"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 4, "N sigma TPC Pion"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 4, "N sigma TPC Proton"};
  Configurable<float> nsigmatpcKa{"nsigmatpcKa", 4, "N sigma TPC Kaon"};
  Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
  Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.05, "DCA bachelor baryon to PV"};
  Configurable<int> mintpccrrows{"mintpccrrows", 50, "min N TPC crossed rows"};
  Configurable<int> dooobrej{"dooobrej", 0, "OOB rejection: 0 no selection, 1 = ITS||TOF, 2 = TOF only for pT > ptthrtof"};
  Configurable<float> ptthrtof{"ptthrtof", 2, "Pt threshold for applying only tof oob rejection"};
  Configurable<float> proplifetime{"proplifetime", 6, "ctau/<ctau>"};
  Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};

  Configurable<bool> doPtDepCosPaCut{"doPtDepCosPaCut", true, "Enable pt dependent cos PA cut"};
  Configurable<bool> doPtDepV0CosPaCut{"doPtDepV0CosPaCut", false, "Enable pt dependent cos PA cut of the V0 daughter"};
  Configurable<bool> doDCAdauToPVCut{"doDCAdauToPVCut", true, "Enable cut DCA daughter track to PV"};
  Configurable<bool> doCascadeCosPaCut{"doCascadeCosPaCut", true, "Enable cos PA cut"};
  Configurable<bool> doV0CosPaCut{"doV0CosPaCut", true, "Enable cos PA cut for the V0 daughter"};
  Configurable<bool> doDCACascadeDauCut{"doDCACascadeDauCut", true, "Enable cut DCA betweenn daughter tracks"};
  Configurable<bool> doDCAV0DauCut{"doDCAV0DauCut", true, "Enable cut DCA betweenn V0 daughter tracks"};
  Configurable<bool> doCascadeRadiusCut{"doCascadeRadiusCut", true, "Enable cut on the cascade radius"};
  Configurable<bool> doV0RadiusCut{"doV0RadiusCut", true, "Enable cut on the V0 radius"};
  Configurable<bool> doDCAV0ToPVCut{"doDCAV0ToPVCut", true, "Enable cut DCA of V0 to PV"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", true, "Enable cut N sigma TPC"};
  Configurable<bool> doCtauCut{"doCtauCut", true, "Enable cut on proplifetime"};
  Configurable<bool> doBachelorBaryonCut{"doBachelorBaryonCut", true, "Enable Bachelor-Baryon cut "};
  Configurable<bool> doProperLifeTimeCut{"doProperLifeTimeCut", true, "Enable proper life-time cut "};

  Partition<soa::Join<aod::CascCores, aod::CascExtras>> negCasc = aod::cascdata::sign < 0;
  Partition<soa::Join<aod::CascCores, aod::CascExtras>> posCasc = aod::cascdata::sign > 0;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101,0,101}});

    histos.add("hCandidate", "hCandidate", HistType::kTH1F, {{22, -0.5, 21.5}});

    TString CutLabel[22] = {"All", "MassWin", "y", "EtaDau", "DCADauToPV", "CascCosPA", "V0CosPA", "DCACascDau", "DCAV0Dau", "rCasc", "rCascMax", "rV0", "rV0Max", "DCAV0ToPV", "LambdaMass", "nSigmaTPCV0Dau", "Bach-baryon", "NTPCrows", "OOBRej", "nSigmaTPCbachelor", "ctau", "CompDecayMass"};
    for (Int_t i = 1; i <= histos.get<TH1>(HIST("hCandidate"))->GetNbinsX(); i++) {
      histos.get<TH1>(HIST("hCandidate"))->GetXaxis()->SetBinLabel(i, CutLabel[i - 1]);
    }

    histos.add("InvMassBefSel/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisXiMass,{101,0,101}});
    histos.add("InvMassBefSel/hPositiveCascade", "hPositiveCascade", {HistType::kTH3F, {axisPt, axisXiMass,{101,0,101}}});

    if(!isXi){
      histos.get<TH3>(HIST("InvMassBefSel/hNegativeCascade"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      histos.get<TH3>(HIST("InvMassBefSel/hPositiveCascade"))->GetYaxis()->Set(200, 1.572f, 1.772f);
    }

    histos.addClone("InvMassBefSel/", "InvMassAfterSel/");
  }

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {

    if (TMath::Abs(coll.posZ()) > zVertexCut) {
      return;
    }
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    if (coll.centFT0C() > 100) {
      return;
    }

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());

    for (auto& casc : Cascades) {

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      // To have trace of how it was before selections
      if (casc.sign() < 0) {
        if (isXi) histos.fill(HIST("InvMassBefSel/hNegativeCascade"), casc.pt(), casc.mXi(), coll.centFT0C());
        else histos.fill(HIST("InvMassBefSel/hNegativeCascade"), casc.pt(), casc.mOmega(), coll.centFT0C());
      }
      if (casc.sign() > 0) {
        if (isXi) histos.fill(HIST("InvMassBefSel/hPositiveCascade"), casc.pt(), casc.mXi(), coll.centFT0C());
        else histos.fill(HIST("InvMassBefSel/hPositiveCascade"), casc.pt(), casc.mOmega(), coll.centFT0C());
      }

      if (isXi) {
        if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) > masswin)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(casc.yXi()) > rapCut)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) > masswin)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
        if (TMath::Abs(casc.yOmega()) > rapCut)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }

      auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});
      if (TMath::Abs(poseta) > etaDauCut || TMath::Abs(negeta) > etaDauCut || TMath::Abs(bacheta) > etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if(doDCAdauToPVCut){
        if (TMath::Abs(casc.dcabachtopv()) < dcaBachToPV)
          continue;
        if(casc.sign() > 0 && (TMath::Abs(casc.dcanegtopv()) < dcaBaryonToPV || TMath::Abs(casc.dcapostopv()) < dcaMesonToPV))
          continue;
        if(casc.sign() < 0 && (TMath::Abs(casc.dcapostopv()) < dcaBaryonToPV || TMath::Abs(casc.dcanegtopv()) < dcaMesonToPV))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if(doCascadeCosPaCut)
      {
        if(doPtDepCosPaCut){
          double ptdepCut = casCosPaPtParameter / casc.pt();
          if (ptdepCut > 0.3 || casc.pt() < 0.5)
            ptdepCut = 0.3;
          if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < TMath::Cos(ptdepCut))
            continue;
        }else if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospa)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if(doV0CosPaCut)
      {
        if(doPtDepV0CosPaCut)
        {
          double ptdepCut = v0CosPaPtParameter / casc.pt();
          if (ptdepCut > 0.3 || casc.pt() < 0.5)
            ptdepCut = 0.3;
          if (casc.casccosPA(casc.x(),casc.y(),casc.z()) < TMath::Cos(ptdepCut))
            continue;
        }else if (casc.casccosPA(casc.x(),casc.y(),casc.z()) < v0cospa)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else ++counter;

      if(doDCACascadeDauCut){
        if (casc.dcacascdaughters() > dcacascdau)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if(doDCAV0DauCut){
        if (casc.dcaV0daughters() > dcav0dau)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if (doCascadeRadiusCut){
        if (casc.cascradius() < minRadius)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);

        if (casc.cascradius() > maxRadius)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else counter+=2;

      if(doV0RadiusCut){
        if (casc.v0radius() < minV0Radius)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
        if (casc.v0radius() > maxV0Radius)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else counter+=2;

      if(doDCAV0ToPVCut){
        if (TMath::Abs(casc.dcav0topv(casc.x(),casc.y(),casc.z())) < dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if (TMath::Abs(casc.mLambda() - pdgDB->Mass(3122)) > lambdaMassWin)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if(doNTPCSigmaCut){
        if (casc.sign() < 0) {
          if( TMath::Abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || TMath::Abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi )
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else if (casc.sign() > 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || TMath::Abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr )
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }
      }else ++counter;

      if(doBachelorBaryonCut){
        if ( (casc.bachBaryonCosPA() > bachBaryonCosPA || TMath::Abs(casc.bachBaryonDCAxyToPV()) < bachBaryonDCAxyToPV)) { // Bach-baryon selection if required
          continue;
        }
        histos.fill(HIST("hCandidate"), ++counter);
      }else ++counter;

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());

      double invmass;
      if (isXi) {
        if(doNTPCSigmaCut){
          if (TMath::Abs(bachExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }else ++counter;

        if(doProperLifeTimeCut){
          if ( pdgDB->Mass(3312) * cascpos / (cascptotmom + 1e-13) > proplifetime * ctauxiPDG)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }else ++counter;

        if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) < rejcomp)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
        invmass = casc.mXi();
      } else {
        if(doNTPCSigmaCut){
          if (TMath::Abs(bachExtra.tpcNSigmaKa()) > nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }else ++counter;

        if(doProperLifeTimeCut){
          if (pdgDB->Mass(3334) * cascpos / (cascptotmom + 1e-13) > proplifetime * ctauomegaPDG)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }else ++counter;
        if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) < rejcomp)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
        invmass = casc.mOmega();
      }

      if(casc.sign() < 0) histos.fill(HIST("InvMassAfterSel/hNegativeCascade"), casc.pt(), invmass, coll.centFT0C());
      else histos.fill(HIST("InvMassAfterSel/hPositiveCascade"), casc.pt(), invmass, coll.centFT0C());

    }
  }

  PROCESS_SWITCH(derivedCascadeAnalysis, processCascades, "cascade analysis, run3 data and rec MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedCascadeAnalysis>(cfgc)};
}
