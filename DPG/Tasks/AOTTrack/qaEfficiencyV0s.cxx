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
/// \file   qaEfficiencyV0s.cxx
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to analyse MC to produce efficiency vs pT, eta and phi of V0s.
///         The efficiency for particles is computed according to the PDG code (sign included and not charge)
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

// ROOT includes
#include "TPDGCode.h"
#include "TEfficiency.h"
#include "THashList.h"

using namespace o2::framework;

struct QaEfficiencyV0s {
  // Particle information
  static constexpr int nSpecies = 1; // One per PDG
  static constexpr const char* particleTitle[nSpecies] = {"K0s"};
  static constexpr int PDGs[nSpecies] = {kK0Short};
  // Particle only selection
  Configurable<bool> doK0s{"do-k0s", false, "Flag to run with the PDG code of k0s"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "Rapidity cut"};
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  OutputObj<THashList> listEfficiencyMC{"EfficiencyMC"};
  // Histograms

  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisPt{ptBins, "#it{p}_{T} GeV/#it{c}"};
    auto h = registry.add<TH1>("Pos/PtNum_310", "Pos/Num_310", kTH1F, {axisPt});
    registry.add("Pos/PtDen_310", "Pos/Den_310", kTH1F, {axisPt});
    registry.add("Neg/PtNum_310", "Neg/Num_310", kTH1F, {axisPt});
    registry.add("Neg/PtDen_310", "Neg/Den_310", kTH1F, {axisPt});
    registry.add("Prm/Pos/PtNum_310", "Prm/Pos/Num_310", kTH1F, {axisPt});
    registry.add("Prm/Pos/PtDen_310", "Prm/Pos/Den_310", kTH1F, {axisPt});
    registry.add("Prm/Neg/PtNum_310", "Prm/Neg/Num_310", kTH1F, {axisPt});
    registry.add("Prm/Neg/PtDen_310", "Prm/Neg/Den_310", kTH1F, {axisPt});
    registry.add("PrmRap/Pos/PtNum_310", "PrmRap/Pos/Num_310", kTH1F, {axisPt});
    registry.add("PrmRap/Pos/PtDen_310", "PrmRap/Pos/Den_310", kTH1F, {axisPt});
    registry.add("PrmRap/Neg/PtNum_310", "PrmRap/Neg/Num_310", kTH1F, {axisPt});
    registry.add("PrmRap/Neg/PtDen_310", "PrmRap/Neg/Den_310", kTH1F, {axisPt});

    TAxis* axis = h->GetXaxis();
    listEfficiencyMC.setObject(new THashList);
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPt_pdg%d", PDGs[0]), Form("efficiencyPt_pdg%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPt_pdg-%d", PDGs[0]), Form("efficiencyPt_pdg-%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrm_pdg%d", PDGs[0]), Form("efficiencyPtPrm_pdg%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrm_pdg-%d", PDGs[0]), Form("efficiencyPtPrm_pdg-%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrmRap_pdg%d", PDGs[0]), Form("efficiencyPtPrmRap_pdg%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrmRap_pdg-%d", PDGs[0]), Form("efficiencyPtPrmRap_pdg-%d", PDGs[0]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
  }

  // MC process
  // Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  // void process(o2::soa::Join<o2::aod::V0Datas, o2::aod::McV0Labels> const& V0s,
  void process(o2::aod::McV0Labels const& V0s,
               o2::aod::McParticles const& mcParticles)
  {
    for (auto const& v0 : V0s) {
      if (v0.has_mcParticle()) {
        auto mcparticle = v0.mcParticle();
        if (mcparticle.pdgCode() == PDGs[0]) {
          registry.fill(HIST("Pos/PtNum_310"), mcparticle.pt());
          if (mcparticle.isPhysicalPrimary()) {
            registry.fill(HIST("Prm/Pos/PtNum_310"), mcparticle.pt());
            if (TMath::Abs(mcparticle.y()) < rapidityCut) {
              registry.fill(HIST("PrmRap/Pos/PtNum_310"), mcparticle.pt());
            }
          }
        } else if (mcparticle.pdgCode() == -PDGs[0]) {
          registry.fill(HIST("Neg/PtNum_310"), mcparticle.pt());
          if (mcparticle.isPhysicalPrimary()) {
            registry.fill(HIST("Prm/Neg/PtNum_310"), mcparticle.pt());
            if (TMath::Abs(mcparticle.y()) < rapidityCut) {
              registry.fill(HIST("PrmRap/Neg/PtNum_310"), mcparticle.pt());
            }
          }
        }
      }
    }
    for (auto const& mcparticle : mcParticles) {
      if (mcparticle.pdgCode() == PDGs[0]) {
        registry.fill(HIST("Pos/PtDen_310"), mcparticle.pt());
        if (mcparticle.isPhysicalPrimary()) {
          registry.fill(HIST("Prm/Pos/PtDen_310"), mcparticle.pt());
          if (TMath::Abs(mcparticle.y()) < rapidityCut) {
            registry.fill(HIST("PrmRap/Pos/PtDen_310"), mcparticle.pt());
          }
        }
      } else if (mcparticle.pdgCode() == -PDGs[0]) {
        registry.fill(HIST("Neg/PtDen_310"), mcparticle.pt());
        if (mcparticle.isPhysicalPrimary()) {
          registry.fill(HIST("Prm/Neg/PtDen_310"), mcparticle.pt());
          if (TMath::Abs(mcparticle.y()) < rapidityCut) {
            registry.fill(HIST("PrmRap/Neg/PtDen_310"), mcparticle.pt());
          }
        }
      }
    }

    static_cast<TEfficiency*>(listEfficiencyMC->At(0))->SetPassedHistogram(*registry.get<TH1>(HIST("Pos/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(0))->SetTotalHistogram(*registry.get<TH1>(HIST("Pos/PtDen_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(1))->SetPassedHistogram(*registry.get<TH1>(HIST("Neg/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(1))->SetTotalHistogram(*registry.get<TH1>(HIST("Neg/PtDen_310")), "f");

    static_cast<TEfficiency*>(listEfficiencyMC->At(2))->SetPassedHistogram(*registry.get<TH1>(HIST("Prm/Pos/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(2))->SetTotalHistogram(*registry.get<TH1>(HIST("Prm/Pos/PtDen_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(3))->SetPassedHistogram(*registry.get<TH1>(HIST("Prm/Neg/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(3))->SetTotalHistogram(*registry.get<TH1>(HIST("Prm/Neg/PtDen_310")), "f");

    static_cast<TEfficiency*>(listEfficiencyMC->At(4))->SetPassedHistogram(*registry.get<TH1>(HIST("PrmRap/Pos/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(4))->SetTotalHistogram(*registry.get<TH1>(HIST("PrmRap/Pos/PtDen_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(5))->SetPassedHistogram(*registry.get<TH1>(HIST("PrmRap/Neg/PtNum_310")), "f");
    static_cast<TEfficiency*>(listEfficiencyMC->At(5))->SetTotalHistogram(*registry.get<TH1>(HIST("PrmRap/Neg/PtDen_310")), "f");
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<QaEfficiencyV0s>(cfgc)}; }
