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
/// \file strangecasctrack.cxx
/// \brief Analysis of strangeness tracking efficiency via primary production of Omega and Xi in Run 3
/// \author Yakiv Paroviak (yakiv.paroviak@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;

struct StrangeCascTrack {

  using TraCascDatas = soa::Join<aod::TraCascIndices, aod::TraCascCores>;
  using CascDatas = soa::Join<aod::CascIndices, aod::CascBBs, aod::CascCores>;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> doProcessPP{"doProcessPP", true, "true for pp, false for PbPb and OO"};
  Configurable<bool> doProcessPbPb{"doProcessPbPb", false, "true for PbPb, false for pp and OO"};
  Configurable<bool> doProcessOO{"doProcessOO", false, "true for OO, false for pp and PbPb"};

  Configurable<bool> doProcessMC{"doProcessMC", false, "true for MC, false for data"};

  Configurable<bool> doRequireFT0{"doRequireFT0", false, "true for offline trigger for Run 3"};
  Configurable<bool> doApplyPID{"doApplyPID", false, "true for offline trigger for Run 3"};
  Configurable<bool> doApplyCuts{"doApplyCuts", true, "true for offline trigger for Run 3"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Multiplicity"};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {2000, 1.6, 1.8}, "#Omega M_{inv} (GeV/c^{2})"};
  ConfigurableAxis axisXiMass{"axisXiMass", {2000, 1.2, 1.4}, "#Xi M_{inv} (GeV/c^{2})"};

  Configurable<double> CutDCAtoPVxy{"CutDCAtoPVxy", 0.02f, "max cascade dca to PV in xy"};
  Configurable<double> CutDCAtoPVz{"CutDCAtoPVz", 0.02f, "max cascade dca to PV in z"};

  void init(InitContext const&)
  {
    histos.add("Events/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});
    histos.add("Events/PVx", "PV x position", kTH1F, {{200, -0.1, 0.1}});
    histos.add("Events/PVy", "PV y position", kTH1F, {{200, -0.1, 0.1}});
    histos.add("Events/PVz", "PV z position", kTH1F, {{100, -20, 20}});
    histos.add("Events/Mult", "Multiplicity", kTH1F, {axisMult});

    histos.add("Tracked/Phi", "Phi", kTH1F, {{100, 0., 2 * M_PI}});
    histos.add("Tracked/Eta", "Eta", kTH1F, {{102, -2.01, 2.01}});
    histos.add("Tracked/DCAxy", "DCA to xy", kTH1F, {{500, 0., 0.5}});
    histos.add("Tracked/DCAz", "DCA to z", kTH1F, {{500, 0., 0.5}});
    histos.add("Tracked/EvMult", "Multiplicity of events with >=1 cascade", kTH1F, {axisMult});
    histos.add("Tracked/MassOmega", "Invariant mass hypothesis", kTH1F, {axisOmegaMass});
    histos.add("Tracked/MassXi", "Invariant mass hypothesis", kTH1F, {axisXiMass});
    histos.add("Tracked/Omega", "", kTHnD, {axisOmegaMass, axisPt, axisMult});
    histos.add("Tracked/Xi", "", kTHnD, {axisXiMass, axisPt, axisMult});

    histos.add("All/Phi", "Phi", kTH1F, {{100, 0., 2 * M_PI}});
    histos.add("All/Eta", "Eta", kTH1F, {{102, -2.01, 2.01}});
    histos.add("All/DCAxy", "DCA to xy", kTH1F, {{1000, 0, 1.}});
    histos.add("All/DCAz", "DCA to z", kTH1F, {{1000, 0, 1.}});
    histos.add("All/EvMult", "Multiplicity of events with >=1 cascade", kTH1F, {axisMult});
    histos.add("All/MassOmega", "Invariant mass hypothesis", kTH1F, {axisOmegaMass});
    histos.add("All/MassXi", "Invariant mass hypothesis", kTH1F, {axisXiMass});
    histos.add("All/Omega", "", kTHnD, {axisOmegaMass, axisPt, axisMult});
    histos.add("All/Xi", "", kTHnD, {axisXiMass, axisPt, axisMult});
  }

  void processTracked(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                      aod::TraCascDatas const& tracascades)
  {
    double mult = doProcessPP ? collision.centFT0M() : collision.centFT0C();
    int64_t casccollid = 0;
    for (auto const& cascade : tracascades) {

      double dcaxy = cascade.dcaXYCascToPV();
      double dcaz = cascade.dcaZCascToPV();
      if (doApplyCuts && ((dcaxy > CutDCAtoPVxy) || (dcaz > CutDCAtoPVz)))
        continue; // DCA check

      if (collision.index() != casccollid) {
        histos.fill(HIST("Tracked/EvMult"), mult); // count and list mult of events with at least one cascade
        casccollid = collision.index();
      }

      double pt = cascade.pt();
      double phi = cascade.phi();
      double eta = cascade.eta();
      double massXi = cascade.mXi();
      double massOmega = cascade.mOmega();

      histos.fill(HIST("Tracked/DCAxy"), dcaxy);
      histos.fill(HIST("Tracked/DCAz"), dcaz);
      histos.fill(HIST("Tracked/Phi"), phi);
      histos.fill(HIST("Tracked/Eta"), eta);
      histos.fill(HIST("Tracked/MassXi"), massXi);
      histos.fill(HIST("Tracked/MassOmega"), massOmega);
      histos.fill(HIST("Tracked/Xi"), massXi, pt, mult);
      histos.fill(HIST("Tracked/Omega"), massOmega, pt, mult);
    }
  }

  void processAll(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                  aod::CascDatas const& cascades)
  {
    histos.fill(HIST("Events/EvCounter"), 0.5);
    double mult = doProcessPP ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("Events/Mult"), mult);
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    histos.fill(HIST("Events/PVx"), pvx);
    histos.fill(HIST("Events/PVy"), pvy);
    histos.fill(HIST("Events/PVz"), pvz);

    int64_t casccollid = 0;
    for (auto const& cascade : cascades) {

      double dcaxy = cascade.dcaXYCascToPV();
      double dcaz = cascade.dcaZCascToPV();
      if (doApplyCuts && ((dcaxy > CutDCAtoPVxy) || (dcaz > CutDCAtoPVz)))
        continue; // DCA check

      if (collision.index() != casccollid) {
        histos.fill(HIST("All/EvMult"), mult); // count and list mult of events with at least one cascade
        casccollid = collision.index();
      }

      double pt = cascade.pt();
      double phi = cascade.phi();
      double eta = cascade.eta();
      double massXi = cascade.mXi();
      double massOmega = cascade.mOmega();

      histos.fill(HIST("All/DCAxy"), dcaxy);
      histos.fill(HIST("All/DCAz"), dcaz);
      histos.fill(HIST("All/Phi"), phi);
      histos.fill(HIST("All/Eta"), eta);
      histos.fill(HIST("All/MassXi"), massXi);
      histos.fill(HIST("All/MassOmega"), massOmega);
      histos.fill(HIST("All/Xi"), massXi, pt, mult);
      histos.fill(HIST("All/Omega"), massOmega, pt, mult);
    }
  }
  PROCESS_SWITCH(StrangeCascTrack, processTracked, "process tracked cascades", true);
  PROCESS_SWITCH(StrangeCascTrack, processAll, "process all cascades", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangeCascTrack>(cfgc),
  };
}
