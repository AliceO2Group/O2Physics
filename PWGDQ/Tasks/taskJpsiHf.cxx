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

/// \file taskJPsiHf.cxx
/// \brief Task for the analysis of J/psi - open HF associate production
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Declarations of various short names
using MyRedEvents = aod::RedJpDmColls;
using MyRedPairCandidatesSelected = aod::RedJpDmDileptons;
using MyRedD0CandidatesSelected = soa::Join<aod::RedJpDmDmesons, aod::RedJpDmD0Masss, aod::RedJpDmDmesBdts>;

struct taskJPsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedDleptDmesAll> redDileptDimesAll;

  // HF configurables
  Configurable<double> massHfCandMin{"massHfCandMin", 1.5, "minimum HF mass"};
  Configurable<double> massHfCandMax{"massHfCandMax", 2.1, "maximum HF mass"};
  // DQ configurables
  Configurable<double> massDileptonCandMin{"massDileptonCandMin", 1, "minimum dilepton mass"};
  Configurable<double> massDileptonCandMax{"massDileptonCandMax", 5, "maximum dilepton mass"};

  void init(o2::framework::InitContext& context)
  {
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <typename TEvent, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TEvent const& event, TDqTrack const& dileptons, THfTrack const& dmesons)
  {
    float ptDilepton = -999;
    float ptDmeson = -999;
    float rapDilepton = -999;
    float rapDmeson = -999;
    float phiDilepton = -999;
    float phiDmeson = -999;
    float deltaRap = -999;
    float deltaPhi = -999;

    for (auto& dilepton : dileptons) {
      ptDilepton = RecoDecay::pt(dilepton.px(), dilepton.py());
      rapDilepton = RecoDecay::y(std::array{dilepton.px(), dilepton.py(), dilepton.pz()}, dilepton.mass());
      phiDilepton = RecoDecay::phi(dilepton.px(), dilepton.py());

      for (auto& dmeson : dmesons) {
        ptDmeson = RecoDecay::pt(dmeson.px(), dmeson.py());
        phiDmeson = RecoDecay::phi(dmeson.px(), dmeson.py());
        deltaPhi = RecoDecay::constrainAngle(phiDilepton - phiDmeson, -o2::constants::math::PIHalf);

        if (dmeson.massD0() > 0) {
          rapDmeson = RecoDecay::y(std::array{dmeson.px(), dmeson.py(), dmeson.pz()}, dmeson.massD0());
          deltaRap = rapDilepton - rapDmeson;
          if ((dilepton.mass() > massDileptonCandMin && dilepton.mass() < massDileptonCandMax) && (dmeson.massD0() > massHfCandMin && dmeson.massD0() < massHfCandMax)) {
            redDileptDimesAll(dilepton.mass(), dmeson.massD0(), ptDilepton, ptDmeson, rapDilepton, rapDmeson, phiDilepton, phiDmeson, deltaRap, deltaPhi, dmeson.bdtBkgMassHypo0(), dmeson.bdtPromptMassHypo0(), dmeson.bdtNonpromptMassHypo0());
          }
        }
        if (dmeson.massD0bar() > 0) {
          rapDmeson = RecoDecay::y(std::array{dmeson.px(), dmeson.py(), dmeson.pz()}, dmeson.massD0bar());
          deltaRap = rapDilepton - rapDmeson;
          if ((dilepton.mass() > massDileptonCandMin && dilepton.mass() < massDileptonCandMax) && (dmeson.massD0() > massHfCandMin && dmeson.massD0() < massHfCandMax)) {
            redDileptDimesAll(dilepton.mass(), dmeson.massD0bar(), ptDilepton, ptDmeson, rapDilepton, rapDmeson, phiDilepton, phiDmeson, deltaRap, deltaPhi, dmeson.bdtBkgMassHypo1(), dmeson.bdtPromptMassHypo1(), dmeson.bdtNonpromptMassHypo1());
          }
        }
      }
    }
  }

  void processRedJspiD0(MyRedEvents::iterator const& event, MyRedPairCandidatesSelected const& dileptons, MyRedD0CandidatesSelected const& dmesons)
  {
    runDileptonDmeson(event, dileptons, dmesons);
  }

  PROCESS_SWITCH(taskJPsiHf, processRedJspiD0, "Process J/psi - D0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskJPsiHf>(cfgc)};
}
