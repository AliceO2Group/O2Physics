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

/// \file taskCorrelationHfeHadrons.cxx
/// \brief HFE-Hadrons azimuthal correlations analysis task - data-like, MC-reco and MC-Gen analyses
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_correlation_electron_hadron;

struct HfTaskCorrelationHfeHadrons {
  // Configurables
  // Deltaphi binning
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 32, "Bins for #Delta#varphi bins"};

  HistogramConfigSpec hCorrelSpec{HistType::kTHnSparseD, {{30, 0., 30.}, {20, 0., 20.}, {nBinsDeltaPhi, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {50, -1.8, 1.8}}};

  HistogramRegistry registry{
    "registry",
    {{"hInclusiveEHCorrel", "Sparse for Delta phi and Delta eta Hadron  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec}}};

  void init(InitContext&)
  {
    registry.get<THnSparse>(HIST("hInclusiveEHCorrel"))->Sumw2();
  }

  // correlation  for electron hadron
  void process(aod::HfEHadronPair const& pairEntries)
  {
    double deltaPhi = -999;
    double deltaEta = -999;
    double ptHadron = -999;
    double ptElectron = -999;

    for (const auto& pairEntry : pairEntries) {

      deltaPhi = pairEntry.deltaPhi();
      deltaEta = pairEntry.deltaEta();
      ptElectron = pairEntry.ptElectron();
      ptHadron = pairEntry.ptHadron();

      registry.fill(HIST("hInclusiveEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationHfeHadrons>(cfgc)};
}
