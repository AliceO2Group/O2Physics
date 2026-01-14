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

// all tbales have the same content (for now)
DECLARE_SOA_TABLE(CJetTFSSs, "AOD", "CJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg); \
DECLARE_SOA_TABLE(CMCDJetTFSSs, "AOD", "CMCDJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg); \
DECLARE_SOA_TABLE(CMCPJetTFSSs, "AOD", "CMCPJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg); \
DECLARE_SOA_TABLE(CEWSJetTFSSs, "AOD", "CEWSJETTFSS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg); \

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETFORMATIONTIMERECLUSTERING_H_



//*********************************************************
//                                                        *
//              Table definitions  old                    *
//                                                        *
//*********************************************************

/*

#ifndef PWGJE_TASKS_JETFORMATIONTIMERECLUSTERING_H_
#define PWGJE_TASKS_JETFORMATIONTIMERECLUSTERING_H_

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h" // IWYU pragma: keep
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"

#include <Framework/ASoA.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace o2::aod
{

namespace jetcollision
{                                                    //!
DECLARE_SOA_COLUMN(PosZ, posZ, float);               //!
DECLARE_SOA_COLUMN(Centrality, centrality, float);   //!
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);     //!
DECLARE_SOA_COLUMN(EventWeight, eventWeight, float); //!
} // namespace jetcollision

namespace jetmccollision
{
DECLARE_SOA_COLUMN(PosZ, posZ, float);               //!
DECLARE_SOA_COLUMN(Accepted, accepted, uint64_t);    //!
DECLARE_SOA_COLUMN(Attempted, attempted, uint64_t);  //!
DECLARE_SOA_COLUMN(XsectGen, xsectGen, float);       //!
DECLARE_SOA_COLUMN(XsectErr, xsectErr, float);       //!
DECLARE_SOA_COLUMN(EventWeight, eventWeight, float); //!
} // namespace jetmccollision

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

namespace splitting
{                                                                                     //!
DECLARE_SOA_COLUMN(Pt, pt, float);                                                    //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                  //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                  //!
DECLARE_SOA_COLUMN(R, r, int);                                                        //!
DECLARE_SOA_COLUMN(SplittingMatchingGeo, splittingMatchingGeo, std::vector<int32_t>); //!
DECLARE_SOA_COLUMN(SplittingMatchingPt, splittingMatchingPt, std::vector<int32_t>);   //!
DECLARE_SOA_COLUMN(SplittingMatchingHF, splittingMatchingHF, std::vector<int32_t>);   //!
} // namespace splitting

// Defines the jet table definition
#define JETSPLITTING_TABLE_DEF(_jet_type_, _jet_description_, _name_, _track_type_, _cand_type_) \
                                                                                                 \
  namespace _name_##splitting                                                                    \
  {                                                                                              \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_##Jet, jet);                                              \
  }                                                                                              \
  namespace _name_##splittingconstituents                                                        \
  {                                                                                              \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(Tracks, tracks, int32_t, _track_type_, "_tracks");       \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(Clusters, clusters, int32_t, JClusters, "_clusters");    \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(Candidates, candidates, int32_t, _cand_type_, "_cand");  \
  }                                                                                              \
  DECLARE_SOA_TABLE(_jet_type_##SPs, "AOD", _jet_description_ "SP",                              \
                    o2::soa::Index<>,                                                            \
                    _name_##splitting::_jet_type_##JetId,                                        \
                    _name_##splittingconstituents::TracksIds,                                    \
                    _name_##splittingconstituents::ClustersIds,                                  \
                    _name_##splittingconstituents::CandidatesIds,                                \
                    splitting::Pt,                                                               \
                    splitting::Eta,                                                              \
                    splitting::Phi,                                                              \
                    splitting::R);

namespace jetoutput
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);                     //!
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);                   //!
DECLARE_SOA_COLUMN(JetEta, jetEta, float);                   //!
DECLARE_SOA_COLUMN(JetY, jetY, float);                       //!
DECLARE_SOA_COLUMN(JetR, jetR, float);                       //!
DECLARE_SOA_COLUMN(JetArea, jetArea, float);                 //!
DECLARE_SOA_COLUMN(JetRho, jetRho, float);                   //!
DECLARE_SOA_COLUMN(JetPerpConeRho, jetPerpConeRho, float);   //!
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, int); //!
} // namespace jetoutput

#define MCCOLL_TABLE_DEF(_jet_type_, _jet_description_, _name_)                                  \
  namespace _name_##mccollisionoutput                                                            \
  {                                                                                              \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, []() -> int { return 0; }); \
  }                                                                                              \
  DECLARE_SOA_TABLE(_jet_type_##MCCOs, "AOD", _jet_description_ "MCCO",                          \
                    jetmccollision::PosZ,                                                        \
                    jetmccollision::Accepted,                                                    \
                    jetmccollision::Attempted,                                                   \
                    jetmccollision::XsectGen,                                                    \
                    jetmccollision::XsectErr,                                                    \
                    jetmccollision::EventWeight,                                                 \
                    _name_##mccollisionoutput::Dummy##_jet_type_<>);

// Defines the jet substrcuture table definition
#define JETSUBSTRUCTURE_TABLE_DEF(_jet_type_, _jet_description_, _name_, _cand_type_, _cand_description_)                                                                                                                                                                                                                                                                                                                                                                     \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  namespace _name_##collisionoutput                                                                                                                                                                                                                                                                                                                                                                                                                                           \
  {                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, []() -> int { return 0; });                                                                                                                                                                                                                                                                                                                                                                              \
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  DECLARE_SOA_TABLE(_jet_type_##COs, "AOD", _jet_description_ "CO", jetcollision::PosZ, jetcollision::Centrality, jetcollision::EventSel, jetcollision::EventWeight, _name_##collisionoutput::Dummy##_jet_type_<>);                                                                                                                                                                                                                                                           \
  using _jet_type_##CO = _jet_type_##COs::iterator;                                                                                                                                                                                                                                                                                                                                                                                                                           \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  namespace _name_##jetoutput                                                                                                                                                                                                                                                                                                                                                                                                                                                 \
  {                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
    DECLARE_SOA_INDEX_COLUMN_CUSTOM(_jet_type_##CO, collision, _jet_description_ "COS");                                                                                                                                                                                                                                                                                                                                                                                      \
    DECLARE_SOA_INDEX_COLUMN_FULL_CUSTOM(Candidate, candidate, int, _cand_type_, _cand_description_ "S", "_0");                                                                                                                                                                                                                                                                                                                                                               \
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
  DECLARE_SOA_TABLE(_jet_type_##Os, "AOD", _jet_description_ "O", _name_##jetoutput::_jet_type_##COId, _name_##jetoutput::CandidateId, jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetoutput::JetY, jetoutput::JetR, jetoutput::JetArea, jetoutput::JetRho, jetoutput::JetPerpConeRho, jetoutput::JetNConstituents);                                                                                                                                              \
  using _jet_type_##O = _jet_type_##Os::iterator;                                                                                                                                                                                                                                                                                                                                                                                                                             \
  namespace _name_##substructure                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  {                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
    DECLARE_SOA_INDEX_COLUMN_CUSTOM(_jet_type_##O, outputTable, _jet_description_ "OS");                                                                                                                                                                                                                                                                                                                                                                                      \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, []() -> int { return 0; });                                                                                                                                                                                                                                                                                                                                                                              \
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
  DECLARE_SOA_TABLE(_jet_type_##SSs, "AOD", _jet_description_ "SS", jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg, _name_##substructure::Dummy##_jet_type_<>); \
  DECLARE_SOA_TABLE(_jet_type_##SSOs, "AOD", _jet_description_ "SSO", _name_##substructure::_jet_type_##OId, jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetTFsubstructure::EnergyMother, jetTFsubstructure::PtLeading, jetTFsubstructure::PtSubLeading, jetTFsubstructure::Theta, jetTFsubstructure::PtLeadingConstituent, jetTFsubstructure::TauForm, jetTFsubstructure::Z, jetTFsubstructure::Ptg, jetTFsubstructure::Thetag, jetTFsubstructure::Zg, jetTFsubstructure::TauFormg);   \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  using _jet_type_##O = _jet_type_##Os::iterator;                                                                                                                                                                                                                                                                                                                                                                                                                             \
  using _jet_type_##SSO = _jet_type_##SSOs::iterator;

#define JETSUBSTRUCTURE_TABLES_DEF(_jet_type_, _jet_description_, _jet_type_full_, _jet_full_description_, _track_type_data_, _cand_type_data_, _cand_description_data_, _track_type_ewsdata_, _cand_type_ewsdata_, _cand_description_ewsdata_, _track_type_mcd_, _cand_type_mcd_, _cand_description_mcd_, _particle_type_, _hfparticle_type_, _hfparticle_description_) \
  JETSUBSTRUCTURE_TABLE_DEF(_jet_type_##Jet, _jet_description_ "JET", _jet_type_##jet, _cand_type_data_, _cand_description_data_)                                                                                                                                                                                                                                        \
  JETSUBSTRUCTURE_TABLE_DEF(_jet_type_##EWSJet, _jet_description_ "EWSJET", _jet_type_##ewsjet, _cand_type_ewsdata_, _cand_description_ewsdata_)                                                                                                                                                                                                                         \
  JETSPLITTING_TABLE_DEF(_jet_type_full_, _jet_description_, _jet_full_description_, _track_type_data_, _cand_type_data_)                                                                                                                                                                                                                                                \
  JETSPLITTING_TABLE_DEF(_jet_type_full_##EventWiseSubtracted, _jet_description_ "EWS", _jet_full_description_##eventwisesubtracted, _cand_type_ewsdata_, _cand_type_ewsdata_)                                                                                                                                                                                           \
  JETSUBSTRUCTURE_TABLE_DEF(_jet_type_##MCDJet, _jet_description_ "MCDJET", _jet_type_##mcdjet, _cand_type_mcd_, _cand_description_mcd_)                                                                                                                                                                                                                                 \
  JETSUBSTRUCTURE_TABLE_DEF(_jet_type_##MCPJet, _jet_description_ "MCPJET", _jet_type_##mcpjet, _hfparticle_type_, _hfparticle_description_)                                                                                                                                                                                                                             \
  MCCOLL_TABLE_DEF(_jet_type_##MCPJet, _jet_description_ "MCPJET", _jet_type_##mcpjet)                                                                                                                                                                                                                                                                                   \
  JETSPLITTING_TABLE_DEF(_jet_type_full_##MCDetectorLevel, _jet_description_ "D", _jet_full_description_##mcdetectorlevel, _track_type_mcd_, _cand_type_mcd_)                                                                                                                                                                                                            \
  JETSPLITTING_TABLE_DEF(_jet_type_full_##MCParticleLevel, _jet_description_ "P", _jet_full_description_##mcparticlelevel, _particle_type_, _hfparticle_type_)

JETSUBSTRUCTURE_TABLES_DEF(C, "C", Charged, charged, JTracks, CJetCOs, "CJETCO", JTrackSubs, CEWSJetCOs, "CEWSJETCO", JTracks, CMCDJetCOs, "CMCDJETCO", JMcParticles, CMCPJetCOs, "CMCPJETCO");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETFORMATIONTIMERECLUSTERING_H_

*/

//*********************************************************
//                                                        *
//              Begin of the task                         *
//                                                        *
//*********************************************************

#include "PWGJE/Core/FastJetUtilities.h"
//#include "PWGJE/Core/JetFinder.h"
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
/*

  Produces<aod::CJetSSs> jetSubstructureDataTable;
  Produces<aod::CMCDJetSSs> jetSubstructureMCDTable;
  Produces<aod::CMCPJetSSs> jetSubstructureMCPTable;
  Produces<aod::CEWSJetSSs> jetSubstructureDataSubTable;
*/
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
      auto theta = parentSubJet1.delta_R(parentSubJet2); // this is deltaR - divide by R in postprocessing
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
    cfgc, TaskName{"jet-formationtimereclustering"})};
}
