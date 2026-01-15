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

// jet analysis tasks (subscribing to jet finder task)
//
/// \author Morgan Knuesel & Johanna Lömker
//

//*********************************************************
//                                                        *
//              Table definitions                         *
//                                                        *
//*********************************************************

#ifndef PWGJE_TASKS_JETFORMATIONTIMERECLUSTERING_H_
#define PWGJE_TASKS_JETFORMATIONTIMERECLUSTERING_H_

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h" // IWYU pragma: keep
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"
#include "PWGJE/DataModel/JetSubstructure.h" // new

#include <Framework/ASoA.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace o2::aod
{
// new part
namespace jetTFsubstructure
{                                                                      //!
DECLARE_SOA_COLUMN(EnergyMother, energyMother, std::vector<float>);    //!
DECLARE_SOA_COLUMN(PtLeading, ptLeading, std::vector<float>);          //!
DECLARE_SOA_COLUMN(PtSubLeading, ptSubLeading, std::vector<float>);    //!
DECLARE_SOA_COLUMN(Theta, theta, std::vector<float>);                  //!
DECLARE_SOA_COLUMN(PtLeadingConstituent, ptLeadingConstituent, float); //!
DECLARE_SOA_COLUMN(TauForm, tauForm, std::vector<float>);              //!

DECLARE_SOA_COLUMN(Z, z, std::vector<float>);               //!
DECLARE_SOA_COLUMN(Ptg, ptg, std::vector<float>);           //!
DECLARE_SOA_COLUMN(Thetag, thetag, std::vector<float>);     //!
DECLARE_SOA_COLUMN(Zg, zg, std::vector<float>);             //!
DECLARE_SOA_COLUMN(TauFormg, tauFormg, std::vector<float>); //!
                                                            //!
} // namespace jetTFsubstructure

// all tables have the same content (for now)
DECLARE_SOA_TABLE(CJetTFSSs, "AOD", "CJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg);
DECLARE_SOA_TABLE(CMCDJetTFSSs, "AOD", "CMCDJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg);
DECLARE_SOA_TABLE(CMCPJetTFSSs, "AOD", "CMCPJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg);
DECLARE_SOA_TABLE(CEWSJetTFSSs, "AOD", "CEWSJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETFORMATIONTIMERECLUSTERING_H_

//*********************************************************
//                                                        *
//              Begin of the task                         *
//                                                        *
//*********************************************************

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetUtilities.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>

#include <cmath>
#include <cstdint>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FormationTimeReclustering {

  Produces<aod::CJetTFSSs> jetSubstructureDataTable;
  Produces<aod::CMCDJetTFSSs> jetSubstructureMCDTable;
  Produces<aod::CMCPJetTFSSs> jetSubstructureMCPTable;
  Produces<aod::CEWSJetTFSSs> jetSubstructureDataSubTable;

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<double> genKTp{"genKTp", 0., "select p value for generalized kT alogrithm"}; // CA: p=0, tau: p=0.5

  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<float> thetaVec;
  std::vector<float> zVec;
  std::vector<float> tauFormVec;

  // groomed vectors:
  std::vector<float> ptgVec;
  std::vector<float> thetagVec;
  std::vector<float> zgVec;
  std::vector<float> taugVec;

  float leadingConstituentPt;
  float ptJet;
  float phiJet;
  float etaJet;
  HistogramRegistry registry;

  void init(InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_tg", ";#it{p}_{T,jet} (GeV/#it{c});#it{#tau}_{g} (fm/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {20, 0.0, 10}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_tg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{#tau}_{g}^{part} (fm/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {20, 0.0, 10}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_tg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{#tau}_{g} (fm/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {20, 0.0, 10}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.fastjetExtraParam = genKTp;                         // in jetfinder we use p = -1 for anti kt jetfinding. Then we do time recl. with p=0.5, kt p =1, ca p=0
    jetReclusterer.algorithm = fastjet::JetAlgorithm::genkt_algorithm; // gen kt is enum 3 in jetfinder setup
  }

  Preslice<aod::JetTracks> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetTracksSub> TracksPerCollisionDataSub = aod::bkgcharged::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  template <bool isMCP, bool isSubtracted, typename T, typename U>
  void jetReclustering(T const& jet, U& splittingTable)
  {
    energyMotherVec.clear();
    ptLeadingVec.clear();
    ptSubLeadingVec.clear();
    thetaVec.clear();
    jetReclustered.clear();
    tauFormVec.clear();
    zVec.clear();
    // groomed
    ptgVec.clear();
    thetagVec.clear();
    zgVec.clear();
    taugVec.clear();

    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;

    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      std::vector<int32_t> tracks;
      std::vector<int32_t> candidates;
      std::vector<int32_t> clusters;
      for (const auto& constituent : sorted_by_pt(parentSubJet2.constituents())) {
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          tracks.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
      }
      splittingTable(jet.globalIndex(), tracks, clusters, candidates, parentSubJet2.perp(), parentSubJet2.eta(), parentSubJet2.phi(), 0);
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);                                                                        // this is deltaR - divide by R in postprocessing
      auto tau = (parentSubJet1.perp() + parentSubJet2.perp()) / (parentSubJet1.perp() * parentSubJet2.perp() * theta * theta); // as in run2 aliphysics
      energyMotherVec.push_back(daughterSubJet.e());
      ptLeadingVec.push_back(parentSubJet1.pt());
      ptSubLeadingVec.push_back(parentSubJet2.pt());
      thetaVec.push_back(theta);
      tauFormVec.push_back(tau);
      zVec.push_back(z);

      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          auto zg = z;
          auto rg = theta;
          auto tg = tau;
          ptgVec.push_back(jet.pt());
          thetagVec.push_back(rg);
          taugVec.push_back(tg);
          zgVec.push_back(zg);
          if constexpr (!isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg);
            registry.fill(HIST("h2_jet_pt_jet_tg"), jet.pt(), tg);
          }
          if constexpr (!isSubtracted && isMCP) {
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg);
            registry.fill(HIST("h2_jet_pt_part_jet_tg_part"), jet.pt(), tg);
          }
          if constexpr (isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"), jet.pt(), rg);
            registry.fill(HIST("h2_jet_pt_jet_tg_eventwiseconstituentsubtracted"), jet.pt(), tg);
          }
          softDropped = true;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
    }
    if constexpr (!isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd);
    }
  }

  template <bool isSubtracted, typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& jet, U const&, V const&, M& outputTable, N& splittingTable)
  {
    jetConstituents.clear();
    ptJet = jet.pt();
    phiJet = jet.phi();
    etaJet = jet.eta();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    jetReclustering<false, isSubtracted>(jet, splittingTable);
    outputTable(ptJet, phiJet, etaJet, energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, leadingConstituentPt, tauFormVec, zVec, ptgVec, thetagVec, zgVec, taugVec);
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(FormationTimeReclustering, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
                              aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable);
  }
  PROCESS_SWITCH(FormationTimeReclustering, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>::iterator const& jet,
                                          aod::JetTracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, TracksPerCollisionDataSub, jetSubstructureDataSubTable, jetSplittingsDataSubTable);
  }
  PROCESS_SWITCH(FormationTimeReclustering, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(typename soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureMCDTable, jetSplittingsMCDTable);
  }
  PROCESS_SWITCH(FormationTimeReclustering, processChargedJetsMCD, "charged jet substructure", false);

  void processChargedJetsMCP(typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet)
  {
    jetConstituents.clear();
    ptJet = jet.pt();
    phiJet = jet.phi();
    etaJet = jet.eta();
    for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    jetReclustering<true, false>(jet, jetSplittingsMCPTable);
    jetSubstructureMCPTable(ptJet, phiJet, etaJet, energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, leadingConstituentPt, tauFormVec, zVec, ptgVec, thetagVec, zgVec, taugVec);
  }
  PROCESS_SWITCH(FormationTimeReclustering, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<FormationTimeReclustering>(
    cfgc, TaskName{"jet-formation-time-reclustering"})};
}
