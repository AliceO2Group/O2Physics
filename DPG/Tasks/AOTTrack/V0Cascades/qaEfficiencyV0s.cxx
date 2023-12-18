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
static constexpr int nSpecies = 2; // One per PDG
// static constexpr const char* particleTitle[nSpecies] = {"K0s", "-K0s"};
static constexpr int PDGs[nSpecies] = {kK0Short, -kK0Short};
int pdgToIndex(int pdg)
{
  for (int i = 0; i < nSpecies; i++) {
    if (pdg == PDGs[i]) {
      return i;
    }
  }
  return -1;
}
static constexpr int kHistoPtNum = 0;
static constexpr int kHistoPtDen = 1;
static constexpr int kHistoTot = 2;

std::shared_ptr<TH1> histograms[nSpecies][kHistoTot] = {{nullptr}};
std::shared_ptr<TH1> histogramsPrm[nSpecies][kHistoTot] = {{nullptr}};
std::shared_ptr<TH1> histogramsPrmRap[nSpecies][kHistoTot] = {{nullptr}};

struct QaEfficiencyV0s {
  // Particle information
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "Rapidity cut"};
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  OutputObj<THashList> listEfficiencyMC{"EfficiencyMC"};
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    const AxisSpec axisPt{ptBins, "#it{p}_{T} GeV/#it{c}"};
    listEfficiencyMC.setObject(new THashList);

    for (int i = 0; i < nSpecies; i++) {
      // MC efficiency (PDG code
      histograms[i][kHistoPtNum] = registry.add<TH1>(Form("Pt/Num_pdg%i", PDGs[i]), Form("Num %i", PDGs[i]), kTH1F, {axisPt});
      histograms[i][kHistoPtDen] = registry.add<TH1>(Form("Pt/Den_pdg%i", PDGs[i]), Form("Den %i", PDGs[i]), kTH1F, {axisPt});

      histogramsPrm[i][kHistoPtNum] = registry.add<TH1>(Form("Pt/Prm/Num_pdg%i", PDGs[i]), Form("Pt Prm Num %i", PDGs[i]), kTH1F, {axisPt});
      histogramsPrm[i][kHistoPtDen] = registry.add<TH1>(Form("Pt/Prm/Den_pdg%i", PDGs[i]), Form("Pt Prm Den %i", PDGs[i]), kTH1F, {axisPt});

      histogramsPrmRap[i][kHistoPtNum] = registry.add<TH1>(Form("Pt/Prm/Rap/PtNum_pdg%i", PDGs[i]), Form("Pt Prm Rap Num %i", PDGs[i]), kTH1F, {axisPt});
      histogramsPrmRap[i][kHistoPtDen] = registry.add<TH1>(Form("Pt/Prm/Rap/PtDen_pdg%i", PDGs[i]), Form("Pt Prm Rap Den %i", PDGs[i]), kTH1F, {axisPt});

      TAxis* axis = histograms[i][0]->GetXaxis();
      listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPt_pdg%d", PDGs[i]), Form("efficiencyPt_pdg%d", PDGs[i]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
      listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrm_pdg%d", PDGs[i]), Form("efficiencyPtPrm_pdg%d", PDGs[i]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
      listEfficiencyMC->Add(new TEfficiency(Form("efficiencyPtPrmRap_pdg%d", PDGs[i]), Form("efficiencyPtPrmRap_pdg%d", PDGs[i]), axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
    }
  }

  // MC process
  // Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  // void process(o2::soa::Join<o2::aod::V0Datas, o2::aod::McV0Labels> const& V0s,
  void process(o2::aod::McV0Labels const& V0s,
               o2::aod::McParticles const& mcParticles)
  {
    for (auto const& v0 : V0s) { // Numerator
      if (v0.has_mcParticle()) {
        auto mcparticle = v0.mcParticle();
        const auto index = pdgToIndex(mcparticle.pdgCode());
        if (index < 0) {
          continue;
        }
        histograms[index][kHistoPtNum]->Fill(mcparticle.pt());
        if (mcparticle.isPhysicalPrimary()) {
          histogramsPrm[index][kHistoPtNum]->Fill(mcparticle.pt());
          if (TMath::Abs(mcparticle.y()) < rapidityCut) {
            histogramsPrmRap[index][kHistoPtNum]->Fill(mcparticle.pt());
          }
        }
      }
    }
    for (auto const& mcparticle : mcParticles) { // Denominator
      const auto index = pdgToIndex(mcparticle.pdgCode());
      if (index < 0) {
        continue;
      }
      histograms[index][kHistoPtDen]->Fill(mcparticle.pt());
      if (mcparticle.isPhysicalPrimary()) {
        histogramsPrm[index][kHistoPtDen]->Fill(mcparticle.pt());
        if (TMath::Abs(mcparticle.y()) < rapidityCut) {
          histogramsPrmRap[index][kHistoPtDen]->Fill(mcparticle.pt());
        }
      }
    }

    for (int i = 0; i < nSpecies; i++) {
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies))->SetPassedHistogram(*histograms[i][kHistoPtNum].get(), "f");
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies))->SetTotalHistogram(*histograms[i][kHistoPtDen].get(), "f");
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies + 1))->SetPassedHistogram(*histogramsPrm[i][kHistoPtNum].get(), "f");
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies + 1))->SetTotalHistogram(*histogramsPrm[i][kHistoPtDen].get(), "f");
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies + 2))->SetPassedHistogram(*histogramsPrmRap[i][kHistoPtNum].get(), "f");
      static_cast<TEfficiency*>(listEfficiencyMC->At(i * nSpecies + 2))->SetTotalHistogram(*histogramsPrmRap[i][kHistoPtDen].get(), "f");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<QaEfficiencyV0s>(cfgc)}; }
