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
/// \brief Step2 of the  LambdaAntiLambdaProducer.cxx
/// \author Akash Raj (akash.raj.john.babu@cern.ch)
/// \author Yash Patley <yash.patley@cern.ch>

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>

#include "CommonConstants/PhysicsConstants.h"

#include <cstdlib>
#include <cmath>
#include <cstdint>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;


void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(
    {"isMC", VariantType::Bool, false,
     {"Run the Monte Carlo analysis workflow"}});
}

#include <Framework/runDataProcessing.h>

//Defining Lambda & Anti-Lambda Table 

namespace o2::aod
{
namespace lambdaorantilambda
{
// Event information
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);

// false = Lambda hypothesis
// true  = anti-Lambda hypothesis
DECLARE_SOA_COLUMN(IsAntiLambdaHypothesis,
                   isAntiLambdaHypothesis,
                   bool);

// V0 kinematics
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);

// Hypothesis-dependent invariant mass
DECLARE_SOA_COLUMN(Mass, mass, float);

// V0 topology
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(CTau, cTau, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);

// Proton for Lambda; antiproton for anti-Lambda
DECLARE_SOA_COLUMN(ProtonPx, protonPx, float);
DECLARE_SOA_COLUMN(ProtonPy, protonPy, float);
DECLARE_SOA_COLUMN(ProtonPz, protonPz, float);
DECLARE_SOA_COLUMN(ProtonPt, protonPt, float);
DECLARE_SOA_COLUMN(ProtonEta, protonEta, float);
DECLARE_SOA_COLUMN(ProtonPhi, protonPhi, float);
DECLARE_SOA_COLUMN(ProtonDcaToPV, protonDcaToPV, float);
DECLARE_SOA_COLUMN(ProtonTrackId, protonTrackId, int64_t);

// Pion information
DECLARE_SOA_COLUMN(PionPx, pionPx, float);
DECLARE_SOA_COLUMN(PionPy, pionPy, float);
DECLARE_SOA_COLUMN(PionPz, pionPz, float);
DECLARE_SOA_COLUMN(PionPt, pionPt, float);
DECLARE_SOA_COLUMN(PionEta, pionEta, float);
DECLARE_SOA_COLUMN(PionPhi, pionPhi, float);
DECLARE_SOA_COLUMN(PionDcaToPV, pionDcaToPV, float);
DECLARE_SOA_COLUMN(PionTrackId, pionTrackId, int64_t);

// Original positive/negative daughter indices
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t);

// Armenteros–Podolanski variables
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(QtArm, qtArm, float);

} // namespace lambdaorantilambda

DECLARE_SOA_TABLE(
  LambdaOrAntiLambdas,
  "AOD",
  "LAMORANTILAM",
  o2::soa::Index<>,

  lambdaorantilambda::CollisionId,
  lambdaorantilambda::Centrality,
  lambdaorantilambda::Multiplicity,
  lambdaorantilambda::Sel8,
  lambdaorantilambda::VertexZ,

  lambdaorantilambda::IsAntiLambdaHypothesis,

  lambdaorantilambda::Px,
  lambdaorantilambda::Py,
  lambdaorantilambda::Pz,
  lambdaorantilambda::Pt,
  lambdaorantilambda::Eta,
  lambdaorantilambda::Phi,
  lambdaorantilambda::Rapidity,

  lambdaorantilambda::Mass,

  lambdaorantilambda::CosPA,
  lambdaorantilambda::CTau,
  lambdaorantilambda::V0Radius,
  lambdaorantilambda::DcaDau,
  lambdaorantilambda::DcaV0ToPV,

  lambdaorantilambda::ProtonPx,
  lambdaorantilambda::ProtonPy,
  lambdaorantilambda::ProtonPz,
  lambdaorantilambda::ProtonPt,
  lambdaorantilambda::ProtonEta,
  lambdaorantilambda::ProtonPhi,
  lambdaorantilambda::ProtonDcaToPV,
  lambdaorantilambda::ProtonTrackId,

  lambdaorantilambda::PionPx,
  lambdaorantilambda::PionPy,
  lambdaorantilambda::PionPz,
  lambdaorantilambda::PionPt,
  lambdaorantilambda::PionEta,
  lambdaorantilambda::PionPhi,
  lambdaorantilambda::PionDcaToPV,
  lambdaorantilambda::PionTrackId,

  lambdaorantilambda::PosTrackId,
  lambdaorantilambda::NegTrackId,

  lambdaorantilambda::Alpha,
  lambdaorantilambda::QtArm);
 

namespace lambdahyperon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);

DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);

DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(Mass, mass, float);

// Proton for Lambda; antiproton for anti-Lambda.
DECLARE_SOA_COLUMN(ProtonPx, protonPx, float);
DECLARE_SOA_COLUMN(ProtonPy, protonPy, float);
DECLARE_SOA_COLUMN(ProtonPz, protonPz, float);

DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t);

DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);
DECLARE_SOA_COLUMN(ProtonDcaToPV, protonDcaToPV, float);
DECLARE_SOA_COLUMN(PionDcaToPV, pionDcaToPV, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(QtArm, qtArm, float);
} // namespace lambdahyperon


DECLARE_SOA_TABLE(Lambdas, "AOD", "LAMBDAS", o2::soa::Index<>,
                  lambdahyperon::CollisionId,
                  lambdahyperon::Centrality, lambdahyperon::Multiplicity,
                  lambdahyperon::Px, lambdahyperon::Py, lambdahyperon::Pz,
                  lambdahyperon::Pt, lambdahyperon::Eta, lambdahyperon::Phi,
                  lambdahyperon::Rapidity, lambdahyperon::Mass,
                  lambdahyperon::ProtonPx, lambdahyperon::ProtonPy, lambdahyperon::ProtonPz,
                  lambdahyperon::PosTrackId, lambdahyperon::NegTrackId,
          lambdahyperon::CosPA,lambdahyperon::DcaDau,lambdahyperon::DcaV0ToPV,
            lambdahyperon::ProtonDcaToPV,lambdahyperon::PionDcaToPV,
            lambdahyperon::V0Radius,lambdahyperon::Alpha,lambdahyperon::QtArm);

DECLARE_SOA_TABLE(AntiLambdas, "AOD", "ANTILAMBDAS", o2::soa::Index<>,
                  lambdahyperon::CollisionId,
                  lambdahyperon::Centrality, lambdahyperon::Multiplicity,
                  lambdahyperon::Px, lambdahyperon::Py, lambdahyperon::Pz,
                  lambdahyperon::Pt, lambdahyperon::Eta, lambdahyperon::Phi,
                  lambdahyperon::Rapidity, lambdahyperon::Mass,
                  lambdahyperon::ProtonPx, lambdahyperon::ProtonPy, lambdahyperon::ProtonPz,
                  lambdahyperon::PosTrackId, lambdahyperon::NegTrackId,
          lambdahyperon::CosPA,lambdahyperon::DcaDau,lambdahyperon::DcaV0ToPV,
            lambdahyperon::ProtonDcaToPV,lambdahyperon::PionDcaToPV,
            lambdahyperon::V0Radius,lambdahyperon::Alpha,lambdahyperon::QtArm);

// =====================================================================
// MC-reconstructed Lambda/anti-Lambda candidates
// One row per reconstructed V0 candidate
// =====================================================================

namespace mcrecolambda
{
// Reconstructed-event information
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);

// MC truth matching
DECLARE_SOA_COLUMN(HasMcMatch, hasMcMatch, bool);
DECLARE_SOA_COLUMN(IsAntiLambda, isAntiLambda, bool);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);

// Generated kinematics of the matched Lambda/anti-Lambda
DECLARE_SOA_COLUMN(GenPx, genPx, float);
DECLARE_SOA_COLUMN(GenPy, genPy, float);
DECLARE_SOA_COLUMN(GenPz, genPz, float);
DECLARE_SOA_COLUMN(GenPt, genPt, float);
DECLARE_SOA_COLUMN(GenEta, genEta, float);
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);
DECLARE_SOA_COLUMN(GenRapidity, genRapidity, float);

// Reconstructed V0 kinematics
DECLARE_SOA_COLUMN(RecoPx, recoPx, float);
DECLARE_SOA_COLUMN(RecoPy, recoPy, float);
DECLARE_SOA_COLUMN(RecoPz, recoPz, float);
DECLARE_SOA_COLUMN(RecoPt, recoPt, float);
DECLARE_SOA_COLUMN(RecoEta, recoEta, float);
DECLARE_SOA_COLUMN(RecoPhi, recoPhi, float);
DECLARE_SOA_COLUMN(RecoRapidity, recoRapidity, float);

//|pTReco-pTGen|
DECLARE_SOA_COLUMN(PtResolution, ptResolution, float);
// Reconstructed invariant-mass hypotheses
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);

// Reconstructed V0 topology
DECLARE_SOA_COLUMN(DcaDaughters, dcaDaughters, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(CTau, cTau, float);

// Proton-Or-AntiProton reconstructed daughter
DECLARE_SOA_COLUMN(ProPx, proPx, float);
DECLARE_SOA_COLUMN(ProPy, proPy, float);
DECLARE_SOA_COLUMN(ProPz, proPz, float);
DECLARE_SOA_COLUMN(ProPt, proPt, float);
DECLARE_SOA_COLUMN(ProEta, proEta, float);
DECLARE_SOA_COLUMN(ProPhi, proPhi, float);
DECLARE_SOA_COLUMN(ProDcaToPV, proDcaToPV, float);
DECLARE_SOA_COLUMN(ProTPCCrossedRows, proTPCCrossedRows, int);
DECLARE_SOA_COLUMN(ProTPCNSigma, proTPCNSigmaPr, float);

// PiePlus-OrPieMinus- reconstructed daughter
DECLARE_SOA_COLUMN(PiePx, piePx, float);
DECLARE_SOA_COLUMN(PiePy, piePy, float);
DECLARE_SOA_COLUMN(PiePz, piePz, float);
DECLARE_SOA_COLUMN(PiePt, piePt, float);
DECLARE_SOA_COLUMN(PieEta, pieEta, float);
DECLARE_SOA_COLUMN(PiePhi, piePhi, float);
DECLARE_SOA_COLUMN(PieDcaToPV, pieDcaToPV, float);
DECLARE_SOA_COLUMN(PieTPCCrossedRows, pieTPCCrossedRows, int);
DECLARE_SOA_COLUMN(PieTPCNSigma, pieTPCNSigmaPr, float);

} // namespace mcrecolambda

DECLARE_SOA_TABLE(
  McRecoLambdaCandidates,
  "AOD",
  "MCRECOLAMBDA",

  o2::soa::Index<>,

  mcrecolambda::CollisionId,
  mcrecolambda::Sel8,
  mcrecolambda::VertexZ,

  mcrecolambda::HasMcMatch,
  mcrecolambda::IsAntiLambda,
  mcrecolambda::IsPhysicalPrimary,

  mcrecolambda::GenPx,
  mcrecolambda::GenPy,
  mcrecolambda::GenPz,
  mcrecolambda::GenPt,
  mcrecolambda::GenEta,
  mcrecolambda::GenPhi,
  mcrecolambda::GenRapidity,

  mcrecolambda::RecoPx,
  mcrecolambda::RecoPy,
  mcrecolambda::RecoPz,
  mcrecolambda::RecoPt,
  mcrecolambda::RecoEta,
  mcrecolambda::RecoPhi,
  mcrecolambda::RecoRapidity,

  mcrecolambda::PtResolution,

  mcrecolambda::MassLambda,

  mcrecolambda::DcaDaughters,
  mcrecolambda::DcaV0ToPV,
  mcrecolambda::CosPA,
  mcrecolambda::V0Radius,
  mcrecolambda::CTau,

  mcrecolambda::ProPx,
  mcrecolambda::ProPy,
  mcrecolambda::ProPz,
  mcrecolambda::ProPt,
  mcrecolambda::ProEta,
  mcrecolambda::ProPhi,
  mcrecolambda::ProDcaToPV,
  mcrecolambda::ProTPCCrossedRows,
  mcrecolambda::ProTPCNSigma,

  mcrecolambda::PiePx,
  mcrecolambda::PiePy,
  mcrecolambda::PiePz,
  mcrecolambda::PiePt,
  mcrecolambda::PieEta,
  mcrecolambda::PiePhi,
  mcrecolambda::PieDcaToPV,
  mcrecolambda::PieTPCCrossedRows,
  mcrecolambda::PieTPCNSigma
);

// =====================================================================
// MC-generated Lambda/anti-Lambda candidates
// One row per generated particle
// =====================================================================

namespace eveselpassmcgenevent
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(NRecoCollisions, nRecoCollisions, int);//Total Number of RecoCollions which is am image of McCollions
DECLARE_SOA_COLUMN(NSel8Collisions, nSel8Collisions, int);//Number of RecoCollions which passed Sel8()
} // namespace mcgenlambda

DECLARE_SOA_TABLE(
  EveSelPassMcGenLambdaCandidates,
  "AOD",
  "MCGENLAMBDA",

  o2::soa::Index<>,

  eveselpassmcgenevent::McCollisionId,
  eveselpassmcgenevent::NRecoCollisions,
  eveselpassmcgenevent::NSel8Collisions);

}// namespace o2::aod

//***********************************************************************************************************
//UnFilterred Lambda-LambdaBar Structure (Only TPCnsigma and TPC min CorssRows fillters are applied )
//***********************************************************************************************************

struct LambdaAntiLambdaProducer {
  
  Produces<aod::LambdaOrAntiLambdas>
    lambdaOrAntiLambdaTable;

  //Produces<aod::Lambdas> lambdaTable;
  //Produces<aod::AntiLambdas> antilambdaTable;

  HistogramRegistry rLambdaOrAntiLambda{"LambdaORAntiLambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //HistogramRegistry rLambda{"Lambdas", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //HistogramRegistry rAntiLambda{"AntiLambdas", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum number of TPC crossed rows"};
Configurable<float> maxTPCNSigma{"maxTPCNSigma", 3.0f, "Maximum absolute TPC PID n-sigma"};

ConfigurableAxis axisInvariantMass{"axisInvariantMass", {100, 1.08f, 1.2f}, "Invariant-mass axis"};
ConfigurableAxis axisLambdaAntiLambdaPt{"axisLambdaPt", {200, 0.f, 10.f}, "Lambda transverse-momentum axis"};
  

void init(InitContext const&){

  // Event information
  rLambdaOrAntiLambda.add("Event/hCentrality", "Centrality;FT0M centrality (%);entries", HistType::kTH1F, {{100, 0.f, 100.f}});
  rLambdaOrAntiLambda.add("Event/hMultiplicity", "Primary-vertex multiplicity;N_{tracks}^{PV};entries", HistType::kTH1F, {{500, 0.f, 500.f}});
  rLambdaOrAntiLambda.add("Candidate/hHypothesis", "Candidate hypothesis;hypothesis;entries", HistType::kTH1F, {{2, -0.5f, 1.5f}});

  // V0 kinematics
  rLambdaOrAntiLambda.add("Candidate/hPx", "V0 p_{x};p_{x} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
  rLambdaOrAntiLambda.add("Candidate/hPy", "V0 p_{y};p_{y} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
  rLambdaOrAntiLambda.add("Candidate/hPz", "V0 p_{z};p_{z} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
  //rLambdaOrAntiLambda.add("Candidate/hPt", "V0 p_{T};p_{T} (GeV/c);entries", HistType::kTH1F, {{200, 0.f, 10.f}});
  rLambdaOrAntiLambda.add("Candidate/hPt", "V0 p_{T};p_{T} (GeV/c);entries", HistType::kTH1F, {axisLambdaAntiLambdaPt});
  rLambdaOrAntiLambda.add("Candidate/hEta", "V0 pseudorapidity;#eta;entries", HistType::kTH1F, {{100, -2.f, 2.f}});
  rLambdaOrAntiLambda.add("Candidate/hPhi", "V0 azimuth;#varphi;entries", HistType::kTH1F, {{72, 0.f, 2.f * M_PI}});
  rLambdaOrAntiLambda.add("Candidate/hRapidity", "V0 rapidity;y;entries", HistType::kTH1F, {{100, -2.f, 2.f}});
  rLambdaOrAntiLambda.add("Candidate/hMass", "Reconstructed invariant mass;m_{p#pi} (GeV/c^{2});entries", HistType::kTH1F, {axisInvariantMass});
  rLambdaOrAntiLambda.add("Candidate/hMassVsPt", "#Lambda invariant mass versus p_{T};" "p_{T} (GeV/c);" "m_{#Lambda} (GeV/c^{2})", HistType::kTH2F, {axisLambdaAntiLambdaPt, axisInvariantMass});
   
  // Topology
  rLambdaOrAntiLambda.add("Topology/hCosPA", "V0 pointing angle;cos(PA);entries", HistType::kTH1F, {{300, 0.97f, 1.f}});
  rLambdaOrAntiLambda.add("Topology/hCTau", "Lambda proper decay length;c#tau (cm);entries", HistType::kTH1F, {{200, 0.f, 100.f}});
  rLambdaOrAntiLambda.add("Topology/hRadius", "V0 transverse decay radius;r_{xy} (cm);entries", HistType::kTH1F, {{200, 0.f, 200.f}});
  rLambdaOrAntiLambda.add("Topology/hDcaDaughters", "DCA between V0 daughters;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 5.f}});
  rLambdaOrAntiLambda.add("Topology/hDcaV0ToPV", "V0 DCA to PV;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 5.f}});

  // Proton or antiproton
  rLambdaOrAntiLambda.add("Proton/hTPCdEdxVsP", "V0 proton/antiproton TPC signal;p_{TPC} (GeV/c);TPC dE/dx (arb. units)", HistType::kTH2F, {{200, 0.f, 10.f}, {500, 0.f, 1000.f}});
  
  rLambdaOrAntiLambda.add("Proton/hPx", "Proton/antiproton p_{x};p_{x} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
  rLambdaOrAntiLambda.add("Proton/hPy", "Proton/antiproton p_{y};p_{y} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
  rLambdaOrAntiLambda.add("Proton/hPz", "Proton/antiproton p_{z};p_{z} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
  rLambdaOrAntiLambda.add("Proton/hPt", "Proton/antiproton p_{T};p_{T} (GeV/c);entries", HistType::kTH1F, {{200, 0.f, 5.f}});
  rLambdaOrAntiLambda.add("Proton/hEta", "Proton/antiproton pseudorapidity;#eta;entries", HistType::kTH1F, {{100, -2.f, 2.f}});
  rLambdaOrAntiLambda.add("Proton/hPhi", "Proton/antiproton azimuth;#varphi;entries", HistType::kTH1F, {{72, 0.f, 2.f * M_PI}});
  rLambdaOrAntiLambda.add("Proton/hDcaToPV", "Proton/antiproton DCA to PV;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 10.f}});

  // Pion
  rLambdaOrAntiLambda.add("Pion/hTPCdEdxVsP", "V0 pion TPC signal;p_{TPC} (GeV/c);TPC dE/dx (arb. units)", HistType::kTH2F, {{200, 0.f, 10.f}, {500, 0.f, 1000.f}});

  rLambdaOrAntiLambda.add("Pion/hPx", "Pion p_{x};p_{x} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
  rLambdaOrAntiLambda.add("Pion/hPy", "Pion p_{y};p_{y} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
  rLambdaOrAntiLambda.add("Pion/hPz", "Pion p_{z};p_{z} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
  rLambdaOrAntiLambda.add("Pion/hPt", "Pion p_{T};p_{T} (GeV/c);entries", HistType::kTH1F, {{200, 0.f, 5.f}});
  rLambdaOrAntiLambda.add("Pion/hEta", "Pion pseudorapidity;#eta;entries", HistType::kTH1F,{{100, -2.f, 2.f}});
  rLambdaOrAntiLambda.add("Pion/hPhi", "Pion azimuth;#varphi;entries", HistType::kTH1F, {{72, 0.f, 2.f * M_PI}});
  rLambdaOrAntiLambda.add("Pion/hDcaToPV", "Pion DCA to PV;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 10.f}});

  // Armenteros-Podolanski
  rLambdaOrAntiLambda.add("Armenteros/hAlpha", "Armenteros #alpha;#alpha;entries", HistType::kTH1F, {{200, -1.f, 1.f}});
  rLambdaOrAntiLambda.add("Armenteros/hQtArm", "Armenteros q_{T};q_{T} (GeV/c);entries", HistType::kTH1F, {{200, 0.f, 0.4f}});
  rLambdaOrAntiLambda.add("Armenteros/hAlphaVsQt", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", HistType::kTH2F, {{200, -1.f, 1.f}, {200, 0.f, 0.4f}});

  auto* hypothesisAxis = rLambdaOrAntiLambda.get<TH1>(HIST("Candidate/hHypothesis"))->GetXaxis();

  hypothesisAxis->SetBinLabel(1, "Lambda");
  hypothesisAxis->SetBinLabel(2, "Anti-Lambda");
  
  }


  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr, aod::pidTPCPi>;
  using CollisionsWithActivity = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;

  void process(aod::V0Datas const& v0s,
               TracksWithPID const&,
               CollisionsWithActivity const&)
    {

      for (auto const& v0 : v0s) {
      auto collision = v0.collision_as<CollisionsWithActivity>();
      auto posTrack = v0.posTrack_as<TracksWithPID>();
      auto negTrack = v0.negTrack_as<TracksWithPID>();

      if (posTrack.tpcNClsCrossedRows() < minTPCCrossedRows || negTrack.tpcNClsCrossedRows() < minTPCCrossedRows) {
      continue;
      }

      const bool isLambda =
        std::abs(posTrack.tpcNSigmaPr()) < maxTPCNSigma &&
        std::abs(negTrack.tpcNSigmaPi()) < maxTPCNSigma;

      const bool isAntiLambda =
        std::abs(negTrack.tpcNSigmaPr()) < maxTPCNSigma &&
        std::abs(posTrack.tpcNSigmaPi()) < maxTPCNSigma;

      // Reject candidates compatible with neither or with both hypotheses.
      if (isLambda == isAntiLambda) {
        continue;
      }

      const auto& protonORantiprotonTrack = isLambda ? posTrack : negTrack;

      const auto& pionTrack = isLambda ? negTrack : posTrack;

      const float mass = isLambda ? v0.mLambda() : v0.mAntiLambda();

      const float protonDcaToPV = isLambda ? std::abs(v0.dcapostopv()) : std::abs(v0.dcanegtopv());

      const float pionDcaToPV = isLambda ? std::abs(v0.dcanegtopv()) : std::abs(v0.dcapostopv());

      const float cTau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;

    

    lambdaOrAntiLambdaTable(
      // Event information
      v0.collisionId(),
      collision.centFT0M(),
      collision.multNTracksPV(),
      collision.sel8(),
      collision.posZ(),

      // false = Lambda; true = anti-Lambda
      isAntiLambda,

      // V0 kinematics
      v0.px(),
      v0.py(),
      v0.pz(),
      v0.pt(),
      v0.eta(),
      v0.phi(),
      v0.yLambda(), // Lambda/anti-Lambda rapidity; v0.y() is the decay position.

      // Reconstructed hypothesis-dependent mass
      mass,

      // V0 topology
      v0.v0cosPA(),
      cTau,
      v0.v0radius(),
      v0.dcaV0daughters(),
      v0.dcav0topv(),

      // Proton for Lambda; antiproton for anti-Lambda
      protonORantiprotonTrack.px(),
      protonORantiprotonTrack.py(),
      protonORantiprotonTrack.pz(),
      protonORantiprotonTrack.pt(),
      protonORantiprotonTrack.eta(),
      protonORantiprotonTrack.phi(),
      protonDcaToPV,
      protonORantiprotonTrack.globalIndex(),

      // Pion
      pionTrack.px(),
      pionTrack.py(),
      pionTrack.pz(),
      pionTrack.pt(),
      pionTrack.eta(),
      pionTrack.phi(),
      pionDcaToPV,
      pionTrack.globalIndex(),

      // Original positive and negative track IDs
      posTrack.globalIndex(),
      negTrack.globalIndex(),

      // Armenteros-Podolanski
      v0.alpha(),
      v0.qtarm() );


  rLambdaOrAntiLambda.fill(HIST("Event/hCentrality"), collision.centFT0M());

  rLambdaOrAntiLambda.fill(HIST("Event/hMultiplicity"), collision.multNTracksPV());

  rLambdaOrAntiLambda.fill(HIST("Candidate/hHypothesis"), isAntiLambda ? 1.f : 0.f);
  rLambdaOrAntiLambda.fill(HIST("Candidate/hPx"), v0.px());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hPy"), v0.py());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hPz"), v0.pz());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hPt"), v0.pt());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hEta"), v0.eta());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hPhi"), v0.phi());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hRapidity"), v0.yLambda());
  rLambdaOrAntiLambda.fill(HIST("Candidate/hMass"), mass);
  rLambdaOrAntiLambda.fill(HIST("Candidate/hMassVsPt"), v0.pt(), mass);

  rLambdaOrAntiLambda.fill(HIST("Topology/hCosPA"), v0.v0cosPA());
  rLambdaOrAntiLambda.fill(HIST("Topology/hCTau"), cTau);
  rLambdaOrAntiLambda.fill(HIST("Topology/hRadius"), v0.v0radius());
  rLambdaOrAntiLambda.fill(HIST("Topology/hDcaDaughters"), v0.dcaV0daughters());
  rLambdaOrAntiLambda.fill(HIST("Topology/hDcaV0ToPV"), v0.dcav0topv());

  rLambdaOrAntiLambda.fill(HIST("Proton/hTPCdEdxVsP"), protonORantiprotonTrack.tpcInnerParam(), protonORantiprotonTrack.tpcSignal());
  rLambdaOrAntiLambda.fill(HIST("Proton/hPx"), protonORantiprotonTrack.px());
  rLambdaOrAntiLambda.fill(HIST("Proton/hPy"), protonORantiprotonTrack.py());
  rLambdaOrAntiLambda.fill(HIST("Proton/hPz"), protonORantiprotonTrack.pz());
  rLambdaOrAntiLambda.fill(HIST("Proton/hPt"), protonORantiprotonTrack.pt());
  rLambdaOrAntiLambda.fill(HIST("Proton/hEta"), protonORantiprotonTrack.eta());
  rLambdaOrAntiLambda.fill(HIST("Proton/hPhi"), protonORantiprotonTrack.phi());
  rLambdaOrAntiLambda.fill(HIST("Proton/hDcaToPV"), protonDcaToPV);

  rLambdaOrAntiLambda.fill(HIST("Pion/hTPCdEdxVsP"), pionTrack.tpcInnerParam(), pionTrack.tpcSignal());
  rLambdaOrAntiLambda.fill(HIST("Pion/hPx"), pionTrack.px());
  rLambdaOrAntiLambda.fill(HIST("Pion/hPy"), pionTrack.py());
  rLambdaOrAntiLambda.fill(HIST("Pion/hPz"), pionTrack.pz());
  rLambdaOrAntiLambda.fill(HIST("Pion/hPt"), pionTrack.pt());
  rLambdaOrAntiLambda.fill(HIST("Pion/hEta"), pionTrack.eta());
  rLambdaOrAntiLambda.fill(HIST("Pion/hPhi"), pionTrack.phi());
  rLambdaOrAntiLambda.fill(HIST("Pion/hDcaToPV"), pionDcaToPV);

  rLambdaOrAntiLambda.fill(HIST("Armenteros/hAlpha"), v0.alpha());
  rLambdaOrAntiLambda.fill(HIST("Armenteros/hQtArm"), v0.qtarm());
  rLambdaOrAntiLambda.fill(HIST("Armenteros/hAlphaVsQt"), v0.alpha(), v0.qtarm());

    }
  }
};

//***********************************************************************************************************
//Filltered Lambda-LambdaBar Structure
//***********************************************************************************************************
struct LambdaAntiLambdaSelector
{
  Produces<aod::Lambdas> lambdaTable;
  Produces<aod::AntiLambdas> antiLambdaTable;

  HistogramRegistry rSelected{"Selected", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};


  ConfigurableAxis axisInvariantMass{"axisInvariantMass", {100, 1.08f, 1.2f}, "Invariant-mass axis"};
  ConfigurableAxis axisPt{"axisPt", {200, 0.f, 10.f}, "transverse-momentum axis"};

  //HistogramRegistry rSelectedLambda{"SelectedLambdas", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //HistogramRegistry rSelectedAntiLambda{"SelectedAntiLambdas", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};


  void addSelectedHistograms(std::string const& prefix)
  {
    // Event information
    rSelected.add((prefix +"/Event/hCentrality").c_str(), "Centrality;FT0M centrality (%);entries", HistType::kTH1F, {{100, 0.f, 100.f}});
    rSelected.add((prefix +"/Event/hMultiplicity").c_str(), "Primary-vertex multiplicity;N_{tracks}^{PV};entries", HistType::kTH1F, {{500, 0.f, 500.f}});

    // Candidate kinematics
    rSelected.add((prefix +"/Candidate/hPx").c_str(), "Candidate p_{x};p_{x} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
    rSelected.add((prefix +"/Candidate/hPy").c_str(), "Candidate p_{y};p_{y} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
    rSelected.add((prefix +"/Candidate/hPz").c_str(), "Candidate p_{z};p_{z} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});
    rSelected.add((prefix +"/Candidate/hPt").c_str(), "Candidate p_{T};p_{T} (GeV/c);entries", HistType::kTH1F,  {axisPt});
    rSelected.add((prefix +"/Candidate/hEta").c_str(), "Candidate pseudorapidity;#eta;entries", HistType::kTH1F, {{100, -2.f, 2.f}});
    rSelected.add((prefix +"/Candidate/hPhi").c_str(), "Candidate azimuth;#varphi;entries", HistType::kTH1F, {{72, 0.f, 6.283185307f}});
    rSelected.add((prefix +"/Candidate/hRapidity").c_str(), "Candidate rapidity;y;entries", HistType::kTH1F, {{100, -1.f, 1.f}});
    rSelected.add((prefix +"/Candidate/hMass").c_str(), "Candidate invariant mass;m_{p#pi} (GeV/c^{2});entries", HistType::kTH1F,  {axisInvariantMass});
    rSelected.add((prefix +"/Candidate/hMassVsPt").c_str(), "#Candidate invariant mass versus p_{T};" "p_{T} (GeV/c);" "m_{#Lambda} (GeV/c^{2})", HistType::kTH2F, {axisPt, axisInvariantMass});
  
    // Proton or antiproton momentum
    rSelected.add((prefix +"/Proton/hPx").c_str(), "Proton/antiproton p_{x};p_{x} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
    rSelected.add((prefix +"/Proton/hPy").c_str(), "Proton/antiproton p_{y};p_{y} (GeV/c);entries", HistType::kTH1F, {{200, -5.f, 5.f}});
    rSelected.add((prefix +"/Proton/hPz").c_str(), "Proton/antiproton p_{z};p_{z} (GeV/c);entries", HistType::kTH1F, {{200, -10.f, 10.f}});

    // Topology
    rSelected.add((prefix +"/Topology/hCosPA").c_str(), "V0 pointing angle;cos(PA);entries", HistType::kTH1F, {{300, 0.97f, 1.f}});
    rSelected.add((prefix +"/Topology/hDcaDaughters").c_str(), "DCA between V0 daughters;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 5.f}});
    rSelected.add((prefix +"/Topology/hDcaV0ToPV").c_str(), "V0 DCA to primary vertex;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 5.f}});
    rSelected.add((prefix +"/Topology/hProtonDcaToPV").c_str(), "Proton/antiproton DCA to PV;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 10.f}});
    rSelected.add((prefix +"/Topology/hPionDcaToPV").c_str(), "Pion DCA to PV;DCA (cm);entries", HistType::kTH1F, {{200, 0.f, 10.f}});
    rSelected.add((prefix +"/Topology/hRadius").c_str(), "V0 transverse decay radius;r_{xy} (cm);entries", HistType::kTH1F, {{200, 0.f, 200.f}});

    // Armenteros-Podolanski
    rSelected.add((prefix +"/Armenteros/hAlpha").c_str(), "Armenteros #alpha;#alpha;entries", HistType::kTH1F, {{200, -1.f, 1.f}});
    rSelected.add((prefix +"/Armenteros/hQtArm").c_str(), "Armenteros q_{T};q_{T} (GeV/c);entries", HistType::kTH1F, {{200, 0.f, 0.4f}});
    rSelected.add((prefix +"/Armenteros/hAlphaVsQt").c_str(), "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", HistType::kTH2F, {{200, -1.f, 1.f}, {200, 0.f, 0.4f}});
  }

  void init(InitContext const&)
  {
    addSelectedHistograms("Lambda");
    addSelectedHistograms("AntiLambda");
  }

  template <typename Candidate>
void fillSelectedHistograms(bool isAntiLambda,
                            Candidate const& candidate)
{
  if (!isAntiLambda) {
    rSelected.fill(HIST("Lambda/Event/hCentrality"), candidate.centrality());
    rSelected.fill(HIST("Lambda/Event/hMultiplicity"), candidate.multiplicity());

    rSelected.fill(HIST("Lambda/Candidate/hPx"), candidate.px());
    rSelected.fill(HIST("Lambda/Candidate/hPy"), candidate.py());
    rSelected.fill(HIST("Lambda/Candidate/hPz"), candidate.pz());
    rSelected.fill(HIST("Lambda/Candidate/hPt"), candidate.pt());
    rSelected.fill(HIST("Lambda/Candidate/hEta"), candidate.eta());
    rSelected.fill(HIST("Lambda/Candidate/hPhi"), candidate.phi());
    rSelected.fill(HIST("Lambda/Candidate/hRapidity"), candidate.rapidity());
    rSelected.fill(HIST("Lambda/Candidate/hMass"), candidate.mass());
    rSelected.fill(HIST("Lambda/Candidate/hMassVsPt"), candidate.pt(), candidate.mass());

    rSelected.fill(HIST("Lambda/Proton/hPx"), candidate.protonPx());
    rSelected.fill(HIST("Lambda/Proton/hPy"), candidate.protonPy());
    rSelected.fill(HIST("Lambda/Proton/hPz"), candidate.protonPz());

    rSelected.fill(HIST("Lambda/Topology/hCosPA"), candidate.cosPA());
    rSelected.fill(HIST("Lambda/Topology/hDcaDaughters"), candidate.dcaDau());
    rSelected.fill(HIST("Lambda/Topology/hDcaV0ToPV"), candidate.dcaV0ToPV());
    rSelected.fill(HIST("Lambda/Topology/hProtonDcaToPV"), candidate.protonDcaToPV());
    rSelected.fill(HIST("Lambda/Topology/hPionDcaToPV"), candidate.pionDcaToPV());
    rSelected.fill(HIST("Lambda/Topology/hRadius"), candidate.v0Radius());

    rSelected.fill(HIST("Lambda/Armenteros/hAlpha"), candidate.alpha());
    rSelected.fill(HIST("Lambda/Armenteros/hQtArm"), candidate.qtArm());
    rSelected.fill(HIST("Lambda/Armenteros/hAlphaVsQt"), candidate.alpha(), candidate.qtArm());
  } else {
    rSelected.fill(HIST("AntiLambda/Event/hCentrality"), candidate.centrality());
    rSelected.fill(HIST("AntiLambda/Event/hMultiplicity"), candidate.multiplicity());

    rSelected.fill(HIST("AntiLambda/Candidate/hPx"), candidate.px());
    rSelected.fill(HIST("AntiLambda/Candidate/hPy"), candidate.py());
    rSelected.fill(HIST("AntiLambda/Candidate/hPz"), candidate.pz());
    rSelected.fill(HIST("AntiLambda/Candidate/hPt"), candidate.pt());
    rSelected.fill(HIST("AntiLambda/Candidate/hEta"), candidate.eta());
    rSelected.fill(HIST("AntiLambda/Candidate/hPhi"), candidate.phi());
    rSelected.fill(HIST("AntiLambda/Candidate/hRapidity"), candidate.rapidity());
    rSelected.fill(HIST("AntiLambda/Candidate/hMass"), candidate.mass());
    rSelected.fill(HIST("AntiLambda/Candidate/hMassVsPt"), candidate.pt(), candidate.mass());

    rSelected.fill(HIST("AntiLambda/Proton/hPx"), candidate.protonPx());
    rSelected.fill(HIST("AntiLambda/Proton/hPy"), candidate.protonPy());
    rSelected.fill(HIST("AntiLambda/Proton/hPz"), candidate.protonPz());

    rSelected.fill(HIST("AntiLambda/Topology/hCosPA"), candidate.cosPA());
    rSelected.fill(HIST("AntiLambda/Topology/hDcaDaughters"), candidate.dcaDau());
    rSelected.fill(HIST("AntiLambda/Topology/hDcaV0ToPV"), candidate.dcaV0ToPV());
    rSelected.fill(HIST("AntiLambda/Topology/hProtonDcaToPV"), candidate.protonDcaToPV());
    rSelected.fill(HIST("AntiLambda/Topology/hPionDcaToPV"), candidate.pionDcaToPV());
    rSelected.fill(HIST("AntiLambda/Topology/hRadius"), candidate.v0Radius());

    rSelected.fill(HIST("AntiLambda/Armenteros/hAlpha"), candidate.alpha());
    rSelected.fill(HIST("AntiLambda/Armenteros/hQtArm"), candidate.qtArm());
    rSelected.fill(HIST("AntiLambda/Armenteros/hAlphaVsQt"), candidate.alpha(), candidate.qtArm());
  }
}
// ================================================================
// Event selection
// ================================================================

  Configurable<float> cutZVertex{ "cutZVertex", 10.f, "Maximum absolute reconstructed vertex z"};

// ================================================================
// Event selection  Filter 
// ================================================================

  Filter selectedLambdaFilter = aod::lambdaorantilambda::sel8 == true &&
                                nabs(aod::lambdaorantilambda::vertexZ) < cutZVertex;

// ================================================================
// Daughter-track selection
// ================================================================

Configurable<float> minDaughterPt{"minDaughterPt", 0.15f, "Minimum daughter-track pT (GeV/c)"};

Configurable<float> maxDaughterEta{"maxDaughterEta", 0.8f, "Maximum absolute daughter-track eta"};



// =====================================================================
// Applying Daughter Filteres 
// =====================================================================


Filter daughterKinematicFilter =
    aod::lambdaorantilambda::protonPt > minDaughterPt &&

    aod::lambdaorantilambda::pionPt > minDaughterPt &&

    nabs(aod::lambdaorantilambda::protonEta) < maxDaughterEta &&

    nabs(aod::lambdaorantilambda::pionEta) < maxDaughterEta;


// ================================================================
// V0 topology selection
// ================================================================


Configurable<float> maxDcaV0Daughters{"maxDcaV0Daughters", 1.0f, "Maximum DCA between V0 daughters (cm)"};

Configurable<float> maxDcaV0ToPV{"maxDcaV0ToPV", 0.1f, "Maximum DCA of V0 to primary vertex (cm)"};

Configurable<float> minDcaProtonToPV{"minDcaProtonToPV", 0.02f, "Minimum proton DCA to primary vertex (cm)"};

Configurable<float> minDcaPionToPV{"minDcaPionToPV", 0.06f, "Minimum pion DCA to primary vertex (cm)"};

Configurable<float> minV0Radius{"minV0Radius", 0.5f, "Minimum transverse V0 decay radius (cm)"};

Configurable<float> maxV0CTau{"maxV0CTau", 30.0f, "Maximum Lambda proper decay length c-tau (cm)"};

Configurable<float> minV0CosPA{ "minV0CosPA", 0.995f, "Minimum cosine of V0 pointing angle"};

// =====================================================================
// Applying V0 topology Filteres 
// =====================================================================

Filter topologyFilter = 
    aod::lambdaorantilambda::dcaDau < maxDcaV0Daughters &&

    aod::lambdaorantilambda::dcaV0ToPV < maxDcaV0ToPV &&

    aod::lambdaorantilambda::protonDcaToPV > minDcaProtonToPV &&

    aod::lambdaorantilambda::pionDcaToPV > minDcaPionToPV &&

    aod::lambdaorantilambda::v0Radius > minV0Radius &&

    aod::lambdaorantilambda::cTau < maxV0CTau &&

    aod::lambdaorantilambda::cosPA > minV0CosPA;

// ================================================================
// Lambda/anti-Lambda invariant-mass selection
// ================================================================


//Configurable<float> minLambdaMass{"minLambdaMass", 1.1081f, "Minimum Lambda invariant mass (GeV/c2)"};

//Configurable<float> maxLambdaMass{"maxLambdaMass", 1.1231f, "Maximum Lambda invariant mass (GeV/c2)"};


// =====================================================================
// Applying invariant-mass Filteres 
// =====================================================================

//Filter lambdaInvariantMassFilter =
//    aod::lambdaorantilambda::mass > minLambdaMass &&
//
//    aod::lambdaorantilambda::mass < maxLambdaMass;


// ================================================================
// Lambda/anti-Lambda kinematic selection
// ================================================================

Configurable<float> minLambdaPt{"minLambdaPt", 0.5f, "Minimum Lambda pT (GeV/c)"};

Configurable<float> maxLambdaPt{"maxLambdaPt", 4.5f, "Maximum Lambda pT (GeV/c)"};

Configurable<float> maxLambdaRapidity{"maxLambdaRapidity", 0.5f, "Maximum absolute Lambda rapidity"};


// =====================================================================
// Applying Lambda/anti-Lambda kinematic Filteres 
// =====================================================================

Filter lambdakinematicFilter =

    aod::lambdaorantilambda::pt > minLambdaPt &&

    aod::lambdaorantilambda::pt < maxLambdaPt &&

    nabs(aod::lambdaorantilambda::rapidity) < maxLambdaRapidity;


  // All compatible filters above are combined with logical AND.
  using SelectedLambdaOrAntiLambdas = soa::Filtered<aod::LambdaOrAntiLambdas>;

  // ==============================================================
  // Lambda and anti-Lambda partitions
  // ==============================================================

  Partition<SelectedLambdaOrAntiLambdas>
    selectedLambdas = aod::lambdaorantilambda::isAntiLambdaHypothesis == false;

  Partition<SelectedLambdaOrAntiLambdas>
    selectedAntiLambdas = aod::lambdaorantilambda::isAntiLambdaHypothesis == true;


void process(SelectedLambdaOrAntiLambdas const&)
  {
    // ============================================================
    // Fill the final Lambda table
    // ============================================================

    for (auto const& candidate : selectedLambdas) {

      fillSelectedHistograms(false, candidate); // Lambda

      lambdaTable(
        candidate.collisionId(),
        candidate.centrality(),
        candidate.multiplicity(),

        candidate.px(),
        candidate.py(),
        candidate.pz(),

        candidate.pt(),
        candidate.eta(),
        candidate.phi(),

        candidate.rapidity(),
        candidate.mass(),

        candidate.protonPx(),
        candidate.protonPy(),
        candidate.protonPz(),

        candidate.posTrackId(),
        candidate.negTrackId(),

        candidate.cosPA(),
        candidate.dcaDau(),
        candidate.dcaV0ToPV(),

        candidate.protonDcaToPV(),
        candidate.pionDcaToPV(),

        candidate.v0Radius(),
        candidate.alpha(),
        candidate.qtArm() );
    }

    // ============================================================
    // Fill the final anti-Lambda table
    // ============================================================

    for (auto const& candidate : selectedAntiLambdas) {

      fillSelectedHistograms(true, candidate);  // anti-Lambda

      antiLambdaTable(
        candidate.collisionId(),
        candidate.centrality(),
        candidate.multiplicity(),

        candidate.px(),
        candidate.py(),
        candidate.pz(),

        candidate.pt(),
        candidate.eta(),
        candidate.phi(),

        candidate.rapidity(),
        candidate.mass(),

        candidate.protonPx(),
        candidate.protonPy(),
        candidate.protonPz(),

        candidate.posTrackId(),
        candidate.negTrackId(),

        candidate.cosPA(),
        candidate.dcaDau(),
        candidate.dcaV0ToPV(),

        candidate.protonDcaToPV(),
        candidate.pionDcaToPV(),

        candidate.v0Radius(),
        candidate.alpha(),
        candidate.qtArm() );

    }
  }


};

//***********************************************************************************************************
// MC Reconstructed and MC Generated Lambda-LambdaBar Table Producer Structure
//***********************************************************************************************************
struct LambdaAntiLambdaMcRecoTableProducer
{
  Produces<aod::McRecoLambdaCandidates> mcRecoTable;

  // ================================================================
  // Daughter-track selection
  // ================================================================

  Configurable<float> maxTPCNSigma{ "maxTPCNSigma", 3.f, "Maximum TPC PID n-sigma"};

  // ================================================================
  // Lambda/anti-Lambda invariant-mass selection
  // ================================================================

  //Configurable<float> minLambdaMass{"minLambdaMass", 1.1081f, "Minimum Lambda invariant mass (GeV/c2)"};
  //Configurable<float> maxLambdaMass{"maxLambdaMass", 1.1231f, "Maximum Lambda invariant mass (GeV/c2)"};

  // =====================================================================
  // Applying invariant-mass Filteres 
  // =====================================================================

  //Filter lambdaInvariantMassFilter =
  //    aod::v0data::mass > minLambdaMass &&
  //    aod::v0data::mass < maxLambdaMass;

  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr, aod::McTrackLabels>;
  using V0CandidatesMC = soa::Join<aod::V0Datas, aod::McV0Labels>;
  //using SelectedV0CandidatesMC = soa::Filtered<V0CandidatesMC>;

  void processMcReco(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    V0CandidatesMC const& V0s, DaughterTracks const&, aod::McParticles const&)
  {

    for (const auto& v0 : V0s) 
    {

        auto posTrack = v0.posTrack_as<DaughterTracks>();
        auto negTrack = v0.negTrack_as<DaughterTracks>();

        const bool isLambda =
          std::abs(posTrack.tpcNSigmaPr()) < maxTPCNSigma &&
          std::abs(negTrack.tpcNSigmaPi()) < maxTPCNSigma;

        const bool isAntiLambda =
          std::abs(negTrack.tpcNSigmaPr()) < maxTPCNSigma &&
          std::abs(posTrack.tpcNSigmaPi()) < maxTPCNSigma;

        if (isLambda == isAntiLambda) {continue;}

      // Hypothesis-dependent daughter assignments
        const auto& protonTrack = isLambda ? posTrack : negTrack;
        const auto& pionTrack = isLambda ? negTrack : posTrack;

        const float protonDcaToPV = isLambda ? std::abs(v0.dcapostopv()) : std::abs(v0.dcanegtopv());
        const float pionDcaToPV = isLambda ? std::abs(v0.dcanegtopv()) : std::abs(v0.dcapostopv());

      // Hypothesis-dependent reconstructed invariant mass
        const float reconstructedMass =
        isLambda ? v0.mLambda() : v0.mAntiLambda();


      // ---------------------------------------------------------------
      // Default MC values for an unmatched reconstructed V0
      // ---------------------------------------------------------------

        const bool hasMcMatch = v0.has_mcParticle();

      // ---------------------------------------------------------------
      // Generated quantities for a truth-matched reconstructed V0
      // ---------------------------------------------------------------

        if (!hasMcMatch) {continue;}
        const auto mcGenParticle = v0.mcParticle();

        const int pdgCode = mcGenParticle.pdgCode();

        // Require reconstructed PID hypothesis to agree with MC truth.
        const bool correctSpecies = (isLambda && pdgCode == 3122) || (isAntiLambda && pdgCode == -3122);
        if (!correctSpecies) {continue;}

        const bool isTrueAntiLambda = ( pdgCode == -3122 );
        const bool isPhysicalPrimary = mcGenParticle.isPhysicalPrimary();

        const float genPx = mcGenParticle.px();
        const float genPy = mcGenParticle.py();
        const float genPz = mcGenParticle.pz();
        const float genPt = mcGenParticle.pt();
        const float genEta = mcGenParticle.eta();
        const float genPhi = mcGenParticle.phi();
        const float genRapidity = mcGenParticle.y();
        

        // ---------------------------------------------------------------
        // Quantities calculated from the reconstructed candidate
        // ---------------------------------------------------------------

        const float cTau = v0.distovertotmom( collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
        const float ptResolution = std::abs( v0.pt() - genPt );
        // ---------------------------------------------------------------
        // Fill exactly one row for every reconstructed V0
        // ---------------------------------------------------------------

        mcRecoTable(
          // Event
          collision.globalIndex(),
          collision.sel8(),
          collision.posZ(),

          // MC truth matching
          hasMcMatch,
          isTrueAntiLambda,
          isPhysicalPrimary,

          // Generated V0 kinematics
          genPx,
          genPy,
          genPz,
          genPt,
          genEta,
          genPhi,
          genRapidity,

          // Reconstructed V0 kinematics
          v0.px(),
          v0.py(),
          v0.pz(),
          v0.pt(),
          v0.eta(),
          v0.phi(),
          v0.yLambda(),


          //|pTReco-pTGen|
          ptResolution,
          // Reconstructed mass hypotheses
          reconstructedMass,

          // Reconstructed topology
          v0.dcaV0daughters(),
          v0.dcav0topv(),
          v0.v0cosPA(),
          v0.v0radius(),
          cTau,

          // Positive daughter
          protonTrack.px(),
          protonTrack.py(),
          protonTrack.pz(),
          protonTrack.pt(),
          protonTrack.eta(),
          protonTrack.phi(),
          protonDcaToPV,
          static_cast<int>(protonTrack.tpcNClsCrossedRows()),
          protonTrack.tpcNSigmaPr(),


          // Negative daughter
          pionTrack.px(),
          pionTrack.py(),
          pionTrack.pz(),
          pionTrack.pt(),
          pionTrack.eta(),
          pionTrack.phi(),
          pionDcaToPV,
          static_cast<int>(pionTrack.tpcNClsCrossedRows()),
          pionTrack.tpcNSigmaPi());
        
      }

    }

Configurable<float> cutZVertexMcGen{"cutZVertexMcGen", 10.f, "Maximum reconstructed vertex z for MCGen event selection"};

using RecoCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

Produces<aod::EveSelPassMcGenLambdaCandidates> mcGenEventSelectedTable;

void processEveSelPassedMcGenTable( aod::McCollisions::iterator const& mcCollision,
                              soa::SmallGroups<RecoCollisionsMC> const& recoCollisions)
{
  int nRecoCollisions = 0;
  int nSel8Collisions = 0;

  for (const auto& recoCollision : recoCollisions) {
    ++nRecoCollisions;

    if (recoCollision.sel8() && std::abs(recoCollision.posZ()) < cutZVertexMcGen) {
      ++nSel8Collisions;
    }
  }

  mcGenEventSelectedTable(
    mcCollision.globalIndex(),
    nRecoCollisions,
    nSel8Collisions);
}

  PROCESS_SWITCH(LambdaAntiLambdaMcRecoTableProducer, processMcReco, "Produce MC reconstructed table", true);
  PROCESS_SWITCH(LambdaAntiLambdaMcRecoTableProducer, processEveSelPassedMcGenTable, "Produce MC generated table", true);

};
//***********************************************************************************************************
// Using Produced Lambda-LambdaBar MC Reconstructed and MC Generated Lambda-LambdaBar For Efficiency Plots
//***********************************************************************************************************
struct LambdaAntiLambdaEfficiencyPlots
{
  HistogramRegistry rEfficiencyPlots{"EfficiencyPlotsy", {},OutputObjHandlingPolicy::AnalysisObject};

  //Defining Plot axis for later purpose
  ConfigurableAxis axisInvariantMass{"axisInvariantMass", {100, 1.08f, 1.2f}, "Invariant-mass axis"};
  ConfigurableAxis axisMcRecoPt{"axisMcRecoPt", {200, 0.f, 10.f}, "McReco transverse-momentum axis"};
  ConfigurableAxis axisGenPt{"axisAssoMcGenPt", {200, 0.f, 10.f}, "McGen transverse-momentum axis"};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {200, 0.f, 2.f}, "|pt_{Reco} - pT_{Gen}|"};
  ConfigurableAxis axisMcRecoEta{"axisEta", {100, -2.f, 2.f}, "McReco eta axis"};
  ConfigurableAxis axisMcGenEta{"axisAssoMcGenEta", {100, -2.f, 2.f}, "McGen eta axis"};


  // ========================Filters======================================
  //~~~~~~~~~~~~~~~~~~~~~~~McRecoFillters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ================================================================
  // Event selection Filters
  // ================================================================

  Configurable<float> cutZVertex{"cutZVertex", 10.f, "Maximum absolute reconstructed vertex z"};

  Filter mcRecoEventFilter = aod::mcrecolambda::sel8 == true &&
                            nabs(aod::mcrecolambda::vertexZ) < cutZVertex &&
                            aod::mcrecolambda::isPhysicalPrimary == true;
  // ================================================================
  // Topological Selection Filters
  // ================================================================

  Configurable<float> maxDcaV0Daughters{"maxDcaV0Daughters", 1.0f, "Maximum DCA between V0 daughters"};
  Configurable<float> maxDcaV0ToPV{"maxDcaV0ToPV", 0.1f, "Maximum V0 DCA to primary vertex"};
  Configurable<float> minDcaProtonToPV{"minDcaProtonToPV", 0.02f, "Minimum proton DCA to primary vertex"};
  Configurable<float> minDcaPionToPV{"minDcaPionToPV", 0.06f, "Minimum pion DCA to primary vertex"};
  Configurable<float> minV0Radius{"minV0Radius", 0.5f, "Minimum V0 decay radius"};
  Configurable<float> maxV0CTau{"maxV0CTau", 30.0f, "Maximum Lambda c-tau"};
  Configurable<float> minV0CosPA{"minV0CosPA", 0.995f, "Minimum V0 cosine of pointing angle"};

  Filter mcRecoTopologyFilter = aod::mcrecolambda::dcaDaughters < maxDcaV0Daughters &&
                                nabs(aod::mcrecolambda::dcaV0ToPV) < maxDcaV0ToPV &&
                                aod::mcrecolambda::proDcaToPV > minDcaProtonToPV &&
                                aod::mcrecolambda::pieDcaToPV > minDcaPionToPV &&
                                aod::mcrecolambda::v0Radius > minV0Radius &&
                                aod::mcrecolambda::cTau < maxV0CTau &&
                                aod::mcrecolambda::cosPA > minV0CosPA;

  // ================================================================
  // kinematics Filter
  // ================================================================

  Configurable<float> minLambdaPt{"minLambdaPt", 0.5f, "Minimum Lambda pT"};
  Configurable<float> maxLambdaPt{"maxLambdaPt", 4.5f, "Maximum Lambda pT"};
  Configurable<float> maxLambdaRapidity{"maxLambdaRapidity", 0.5f, "Maximum absolute Lambda rapidity"};

  Filter mcRecoKinematicFilter = aod::mcrecolambda::recoPt > minLambdaPt &&
                                aod::mcrecolambda::recoPt < maxLambdaPt &&
                                nabs(aod::mcrecolambda::recoRapidity) < maxLambdaRapidity;

  // ================================================================
  // Daughter-track selection
  // ================================================================

  Configurable<float> minDaughterPt{"minDaughterPt", 0.15f, "Minimum daughter-track pT"};
  Configurable<float> maxDaughterEta{"maxDaughterEta", 0.8f, "Maximum absolute daughter-track eta"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum TPC crossed rows"};

  Filter mcRecoDaughterFilter = aod::mcrecolambda::proPt > minDaughterPt 
                                && aod::mcrecolambda::piePt > minDaughterPt &&
                                nabs(aod::mcrecolambda::proEta) < maxDaughterEta &&
                                nabs(aod::mcrecolambda::pieEta) < maxDaughterEta &&
                                aod::mcrecolambda::proTPCCrossedRows >= minTPCCrossedRows &&
                                aod::mcrecolambda::pieTPCCrossedRows >= minTPCCrossedRows;


  void addEfficiencyPlotsHistograms(std::string const& particle)
  {
    const std::string recoDirectory = particle + "/Reco/";
    const std::string genDirectory = particle + "/Gen/";

    rEfficiencyPlots.add( (recoDirectory + "hMcRecoPt").c_str(), "Reconstructed p_{T};p_{T} (GeV/c);Entries", HistType::kTH1F, {axisMcRecoPt});
    rEfficiencyPlots.add( (recoDirectory + "hAssoMcGenPt").c_str(), "Reco Accociated Gen p_{T};p_{T} (GeV/c);Entries", HistType::kTH1F, {axisGenPt});
    rEfficiencyPlots.add( (recoDirectory + "hMcRecoEta").c_str(), "Reconstructed #eta;#eta;Entries", HistType::kTH1F, {axisMcRecoEta});
    rEfficiencyPlots.add( (recoDirectory + "hAssoMcGenEta").c_str(), "Reco Accociated Gen #eta;#eta;Entries", HistType::kTH1F, {axisMcGenEta});
    rEfficiencyPlots.add( (recoDirectory + "hPtResolution").c_str(), "|pT_{Reco} - pT_{Gen}|", HistType::kTH1F, {axisDeltaPt});
    rEfficiencyPlots.add( (recoDirectory + "hInvariantMass").c_str(), "Reconstructed invariant mass;m_{p#pi} (GeV/c^{2});entries", HistType::kTH1F, {axisInvariantMass});

    rEfficiencyPlots.add( (genDirectory + "hPt").c_str(), "Generated p_{T};p_{T} (GeV/c);Entries", HistType::kTH1F, {axisGenPt});
    rEfficiencyPlots.add( (genDirectory + "hEta").c_str(), "Generated #eta;#eta;Entries", HistType::kTH1F, {axisMcGenEta});
    rEfficiencyPlots.add( (genDirectory + "PrimarAndNonPrimary").c_str(), "Generated Primary And Non-Primary Genp_{T};Genp_{T} (GeV/c);Entries", HistType::kTH1F, {axisGenPt});
    
    
  }

  void init(InitContext const&){
  addEfficiencyPlotsHistograms("Lambda");
  addEfficiencyPlotsHistograms("AntiLambda");
  }

  using SelectedMcRecoLambdas = soa::Filtered<aod::McRecoLambdaCandidates>;

  Partition<SelectedMcRecoLambdas> selectedMcRecoLambdas = aod::mcrecolambda::isAntiLambda == false;
  Partition<SelectedMcRecoLambdas> selectedMcRecoAntiLambdas = aod::mcrecolambda::isAntiLambda == true;

  void processEfficiencyMcReco( SelectedMcRecoLambdas const&)
  {
    for (const auto& lambda : selectedMcRecoLambdas) {
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hMcRecoPt"), lambda.recoPt());
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hAssoMcGenPt"), lambda.genPt());
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hMcRecoEta"), lambda.recoEta());
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hAssoMcGenEta"), lambda.genEta());
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hPtResolution"), lambda.ptResolution());
      rEfficiencyPlots.fill(HIST("Lambda/Reco/hInvariantMass"), lambda.massLambda());
      }

    for (const auto& antiLambda : selectedMcRecoAntiLambdas) {
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hMcRecoPt"), antiLambda.recoPt());
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hAssoMcGenPt"), antiLambda.genPt());
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hMcRecoEta"), antiLambda.recoEta());
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hAssoMcGenEta"), antiLambda.genEta());
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hPtResolution"), antiLambda.ptResolution());
      rEfficiencyPlots.fill(HIST("AntiLambda/Reco/hInvariantMass"), antiLambda.massLambda());
      }
}
  PROCESS_SWITCH(LambdaAntiLambdaEfficiencyPlots, processEfficiencyMcReco, "Fill MC reconstructed efficiency histograms", true);
  
  //McGen Filters
  Filter mcGenEventFilter = aod::eveselpassmcgenevent::nSel8Collisions > 0;

  using SelectedMcGenEvents = soa::Filtered<aod::EveSelPassMcGenLambdaCandidates>;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  void processEfficiencyMcGen(SelectedMcGenEvents::iterator const& selectedMcEvent, aod::McParticles const& allMcParticles)
  {
  // We will slice the particles here later.
    const auto particlesFromThisMcCollision = allMcParticles.sliceBy(mcParticlesPerMcCollision, selectedMcEvent.mcCollisionId());

  // Reconstructed-collision-opportunity weight
    const float weight = static_cast<float>(selectedMcEvent.nSel8Collisions());

    for (const auto& particle : particlesFromThisMcCollision) {

     // Lambda or anti-Lambda
      //if (std::abs(particle.pdgCode()) != 3122) {continue;}
      
      // Generated kinematic acceptance
      if (particle.pt() <= minLambdaPt ||particle.pt() >= maxLambdaPt) {continue;}
      if (std::abs(particle.y()) >= maxLambdaRapidity) {continue;}

      if (particle.pdgCode() == 3122) {
        // Inclusive: primary + non-primary Lambda
        rEfficiencyPlots.fill(HIST("Lambda/Gen/PrimarAndNonPrimary"), particle.pt(), weight);

        if (!particle.isPhysicalPrimary()) {continue;}

        rEfficiencyPlots.fill(HIST("Lambda/Gen/hPt"), particle.pt(), weight);
        rEfficiencyPlots.fill(HIST("Lambda/Gen/hEta"), particle.eta(), weight);

      }else if (particle.pdgCode() == -3122) {
        // Inclusive: primary + non-primary anti-Lambda
        rEfficiencyPlots.fill(HIST("AntiLambda/Gen/PrimarAndNonPrimary"), particle.pt(), weight);

        if (!particle.isPhysicalPrimary()) {continue;}

        rEfficiencyPlots.fill(HIST("AntiLambda/Gen/hPt"), particle.pt(), weight);
        rEfficiencyPlots.fill(HIST("AntiLambda/Gen/hEta"), particle.eta(), weight);
      }
      
    }
  }
  PROCESS_SWITCH(LambdaAntiLambdaEfficiencyPlots, processEfficiencyMcGen, "Fill MC-generated efficiency histograms", true);
};

struct LambdaAntiLambdaSelectionCutFlow
{
  HistogramRegistry cutFlowRegistry{"cutFlowRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event selection
  Configurable<float> cutZVertex{"cutFlowCutZVertex", 10.f, "Maximum absolute reconstructed vertex z"};

  // Daughter selection
  Configurable<float> minDaughterPt{"cutFlowMinDaughterPt", 0.15f, "Minimum daughter pT"};
  Configurable<float> maxDaughterEta{"cutFlowMaxDaughterEta", 0.8f, "Maximum daughter absolute eta"};
  Configurable<int> minTPCCrossedRows{"cutFlowMinTPCCrossedRows", 70, "Minimum TPC crossed rows"};

  // Topology
  Configurable<float> maxDcaV0Daughters{"cutFlowMaxDcaV0Daughters", 1.f, "Maximum daughter DCA"};
  Configurable<float> maxDcaV0ToPV{"cutFlowMaxDcaV0ToPV", 0.1f, "Maximum V0 DCA to PV"};
  Configurable<float> minDcaProtonToPV{"cutFlowMinDcaProtonToPV", 0.02f, "Minimum proton DCA to PV"};
  Configurable<float> minDcaPionToPV{"cutFlowMinDcaPionToPV", 0.06f, "Minimum pion DCA to PV"};
  Configurable<float> minV0Radius{"cutFlowMinV0Radius", 0.5f, "Minimum V0 radius"};
  Configurable<float> maxV0CTau{"cutFlowMaxV0CTau", 30.f, "Maximum Lambda c-tau"};
  Configurable<float> minV0CosPA{"cutFlowMinV0CosPA", 0.995f, "Minimum cosine of pointing angle"};

  // Lambda kinematics
  Configurable<float> minLambdaPt{"cutFlowMinLambdaPt", 0.5f, "Minimum Lambda pT"};
  Configurable<float> maxLambdaPt{"cutFlowMaxLambdaPt", 4.5f, "Maximum Lambda pT"};
  Configurable<float> maxLambdaRapidity{"cutFlowMaxLambdaRapidity", 0.5f, "Maximum absolute Lambda rapidity"};


    void init(InitContext const&)
  {
    AxisSpec eventCutAxis{ 3, -0.5, 2.5, "Event-selection stage"};

    AxisSpec lambdaCutAxis{ 8, -0.5, 7.5, "Candidate-selection stage"};

    cutFlowRegistry.add("Event/hCutFlow", "Event cut flow;Selection;Collisions", HistType::kTH1D, {eventCutAxis});
    cutFlowRegistry.add("Lambda/hCutFlow", "Lambda cut flow;Selection;Candidates", HistType::kTH1D, {lambdaCutAxis});
    cutFlowRegistry.add("AntiLambda/hCutFlow", "Anti-Lambda cut flow;Selection;Candidates", HistType::kTH1D, {lambdaCutAxis});

    auto eventHistogram = cutFlowRegistry.get<TH1>(HIST("Event/hCutFlow"));

    eventHistogram->GetXaxis()->SetBinLabel(1, "All collisions");
    eventHistogram->GetXaxis()->SetBinLabel(2, "sel8");
    eventHistogram->GetXaxis()->SetBinLabel(3, "Vertex z");

    auto lambdaHistogram = cutFlowRegistry.get<TH1>(HIST("Lambda/hCutFlow"));

    auto antiLambdaHistogram =cutFlowRegistry.get<TH1>(HIST("AntiLambda/hCutFlow"));

    for (const auto& histogram : {lambdaHistogram, antiLambdaHistogram}) {
      histogram->GetXaxis()->SetBinLabel(1, "Truth matched + PID selected");
      histogram->GetXaxis()->SetBinLabel(2, "sel8");
      histogram->GetXaxis()->SetBinLabel(3, "Vertex z");
      histogram->GetXaxis()->SetBinLabel(4, "Physical primary");
      histogram->GetXaxis()->SetBinLabel(5, "Daughter acceptance");
      histogram->GetXaxis()->SetBinLabel(6, "TPC crossed rows");
      histogram->GetXaxis()->SetBinLabel(7, "V0 topology");
      histogram->GetXaxis()->SetBinLabel(8, "Lambda kinematics");
    }
  }

    using ReconstructedEvents = soa::Join<aod::Collisions, aod::EvSels>;

  void processEventCutFlow(ReconstructedEvents::iterator const& collision)
  {
    cutFlowRegistry.fill( HIST("Event/hCutFlow"), 0.f);

    if (!collision.sel8()) {return;}

    cutFlowRegistry.fill(HIST("Event/hCutFlow"), 1.f);

    if (std::abs(collision.posZ()) >= cutZVertex) {return;}

    cutFlowRegistry.fill(HIST("Event/hCutFlow"), 2.f);
  }

    void processLambdaCutFlow(aod::McRecoLambdaCandidates const& candidates)
  {
    for (const auto& candidate : candidates) {
      const bool isAntiLambda = candidate.isAntiLambda();

      auto fillStage = [&](float stage) {
        if (isAntiLambda) {
          cutFlowRegistry.fill(HIST("AntiLambda/hCutFlow"), stage);
        } else {
          cutFlowRegistry.fill(HIST("Lambda/hCutFlow"), stage);
        }
      };

      // PID, truth match and correct species were already
      // required by LambdaAntiLambdaMcRecoTableProducer.
      fillStage(0.f);

      if (!candidate.sel8()) {continue;}
      fillStage(1.f);

      if (std::abs(candidate.vertexZ()) >= cutZVertex) {continue;}
      fillStage(2.f);

      if (!candidate.isPhysicalPrimary()) {continue;}
      fillStage(3.f);

      if (candidate.proPt() <= minDaughterPt ||candidate.piePt() <= minDaughterPt ||
          std::abs(candidate.proEta()) >= maxDaughterEta ||
          std::abs(candidate.pieEta()) >= maxDaughterEta) {continue;}
      fillStage(4.f);

      if (candidate.proTPCCrossedRows() < minTPCCrossedRows ||
          candidate.pieTPCCrossedRows() < minTPCCrossedRows) {continue;}
      fillStage(5.f);

      if (candidate.dcaDaughters() >= maxDcaV0Daughters || std::abs(candidate.dcaV0ToPV()) >= maxDcaV0ToPV ||
          candidate.proDcaToPV() <= minDcaProtonToPV || candidate.pieDcaToPV() <= minDcaPionToPV ||
          candidate.v0Radius() <= minV0Radius || candidate.cTau() >= maxV0CTau ||
          candidate.cosPA() <= minV0CosPA) {continue;}
      fillStage(6.f);

      if (candidate.recoPt() <= minLambdaPt || candidate.recoPt() >= maxLambdaPt ||
          std::abs(candidate.recoRapidity()) >= maxLambdaRapidity) {continue;}
      fillStage(7.f);
    }
  }

  PROCESS_SWITCH(LambdaAntiLambdaSelectionCutFlow, processEventCutFlow, "Fill reconstructed-event cut flow", true);
  PROCESS_SWITCH(LambdaAntiLambdaSelectionCutFlow, processLambdaCutFlow, "Fill MCReco Lambda cut flow", true);
};

struct LambdaAntiLambdaPairAnalysis
{
  //Defining Histogram registry
  HistogramRegistry rSpinAnalysis{"SpinPairAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //defining Axis for histograms
  ConfigurableAxis axisDeltaEta{ "axisDeltaEta", {100, -2.f, 2.f}, "#Delta#eta axis"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -M_PI, M_PI}, "#Delta#phi axis"};
  ConfigurableAxis axisDeltaPt{ "axisDeltaPt", {400, -10.f, 10.f}, "#Delta#p_{T} axis"};

  ConfigurableAxis axisPt{"axisPt", {200, 0.f, 10.f}, "transverse-momentum axis"};
  ConfigurableAxis axisEta{"axisEta", {100, -2.f, 2.f}, "#eta axis"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0, 2*M_PI}, "#phi axis"};

  ConfigurableAxis axisInvariantMass{"axisInvariantMass", {100, 1.08f, 1.2f}, "Invariant-mass axis"};
  ConfigurableAxis axisDeltaMass{"axisDeltaInvariantMass", {100, -0.12f, 0.12f}, "Invariant-mass axis"};
  //ConfigurableAxis axisNumberOfPairs{"axisNumberOfPairs", {1001, -0.5f, 1000.5f}, "NumberOfPairs axis"};

  ConfigurableAxis axisCos{"axisCos", {200, -1.f, 1.f}, "#Cos#theta axis"};
  ConfigurableAxis axisMult{"axisMult", {200, 0.f, 200.f}, "Multiplicity axis"};
  ConfigurableAxis axisDeltaMult{"axisDeltaMult", {400, -200.f, 200.f}, "Delta Multiplicity axis"};
  ConfigurableAxis axisCent{"axisCent", {100, 0.f, 100.f}, "Centrality axis"};
  ConfigurableAxis axisDeltaCent{"axisDeltaCent", {200, -100.f, 100.f}, "Delta Centrality axis"};

  //Mixed-Event Compatablity variables
  Configurable<float> compatibilityDeltaPt{"compatibilityDeltaPt", 0.1f, "compatibility pT difference for mixed-candidate matching"};
  Configurable<float> compatibilityDeltaPhi{"compatibilityDeltaPhi", 0.1f, "compatibility phi difference for mixed-candidate matching"};
  Configurable<float> compatibilityDeltaRapidity{"compatibilityDeltaRapidity", 0.1f, "compatibility rapidity difference for mixed-candidate matching"};
  

    void addHistograms(TString const& candidate1, TString const& candidate2){
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Same-Event~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //EventConter histogram
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/EventCount",  "Same-event" + candidate1 + candidate2 + "pairs;" "SameEvent", HistType::kTH1F, {{1,-0.5, 0.5}});
      rSpinAnalysis.add(candidate1 + candidate2 + "/hSameEvent/NumberOfPairs", "Same-event pair count;Counter;Accepted pairs", HistType::kTH1F, {{1,-0.5, 0.5}});
      //Delta Kinaematics histos
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaEtaDelatPhi",  "Same-event" + candidate1 + candidate2 + "pairs;" "#Delta#eta;" "#Delta#phi", HistType::kTH2F, { axisDeltaEta, axisDeltaPhi });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaPtDelatEta",  "Same-event" + candidate1 + candidate2 + "pairs;" "#Delta#p_{T};" "#Delta#eta", HistType::kTH2F, { axisDeltaPt, axisDeltaEta });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaPtDelatPhi",  "Same-event" + candidate1 + candidate2 + "pairs;" "#Delta#p_{T};" "#Delta#phi", HistType::kTH2F, { axisDeltaPt, axisDeltaPhi });
      //Invarient mass and Cent histos
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaMassDeltaCent",  "Same-event" + candidate1 + candidate2 + "pairs;" "#Delta m_{#pi#p} ;" "DeltaCent", HistType::kTH2F, { axisDeltaMass, axisDeltaCent });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaMassDeltaPt",  "Same-event" + candidate1 + candidate2 + "pairs;" "#Delta m_{#pi#p} ;" "Deltap_{T}", HistType::kTH2F, { axisDeltaMass, axisDeltaPt });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/DeltaCentDeltaMult",  "Same-event" + candidate1 + candidate2 + "pairs;" "DeltaCent ;" "DeltaMult", HistType::kTH2F, { axisDeltaCent, axisDeltaMult });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/Mass",  "Same-event" + candidate1 + candidate2 + "pairs;" "Mass-Lambda ;" "Mass-AntiLambda", HistType::kTH2F, { axisInvariantMass, axisInvariantMass });
      
      //Spin histos
        //2d
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/CosMult",  "Same-event" + candidate1 + candidate2 + "pairs;" "cos#theta;" "Mult", HistType::kTH2F, { axisCos, axisMult });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/CosCent",  "Same-event" + candidate1 + candidate2 + "pairs;" "cos#theta;" "Cent", HistType::kTH2F, { axisCos, axisCent });
        //1d
      rSpinAnalysis.add(candidate1 + candidate2+ "/hSameEvent/Cos",  "Same-event" + candidate1 + candidate2 + "pairs;" "cos#theta" , HistType::kTH1F, { axisCos });
      

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Mixed-Event~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //EventConter histogram
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/EventCount",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "MixedEvent", HistType::kTH1F, {{1,-0.5, 0.5}});
      rSpinAnalysis.add(candidate1 + candidate2 + "/hMixedEvent/NumberOfPairs", "Mixed-event pair count;Counter;Accepted pairs", HistType::kTH1F, {{1,-0.5, 0.5}});
      
      //Delta Kinaematics histos
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaEtaDelatPhi",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "#Delta#eta;" "#Delta#phi", HistType::kTH2F, { axisDeltaEta, axisDeltaPhi });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaPtDelatEta",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "#Delta#p_{T};" "#Delta#eta", HistType::kTH2F, { axisDeltaPt, axisDeltaEta });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaPtDelatPhi",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "#Delta#p_{T};" "#Delta#phi", HistType::kTH2F, { axisDeltaPt, axisDeltaPhi });
      //Invarient mass and Cent histos
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaMassDeltaCent",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "#Delta m_{#pi#p} ;" "DeltaCent", HistType::kTH2F, { axisDeltaMass, axisDeltaCent });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaMassDeltaPt",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "#Delta m_{#pi#p} ;" "Deltap_{T}", HistType::kTH2F, { axisDeltaMass, axisDeltaPt });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/DeltaCentDeltaMult",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "DeltaCent ;" "DeltaMult", HistType::kTH2F, { axisDeltaCent, axisDeltaMult });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/Mass",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "Mass-Lambda ;" "Mass-AntiLambda", HistType::kTH2F, { axisInvariantMass, axisInvariantMass });
      
      //Spin histos
        //2d
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/CosMult",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "cos#theta;" "Mult"+candidate1+ ";" "Mult"+candidate2 , HistType::kTH3F, { axisCos, axisMult, axisMult });
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/CosCent",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "cos#theta;" "Cent"+candidate1+ ";" "Cent"+candidate2, HistType::kTH3F, { axisCos, axisCent, axisCent });
        //1d
      rSpinAnalysis.add(candidate1 + candidate2+ "/hMixedEvent/Cos",  "Mixed-event" + candidate1 + candidate2 + "pairs;" "cos#theta" , HistType::kTH1F, { axisCos });
      
    }

  void init(InitContext const&)
    { 
      addHistograms("Lambda", "AntiLambda");
      addHistograms("Lambda", "Lambda");
      addHistograms("AntiLambda", "AntiLambda");
    }

    ROOT::Math::PxPyPzMVector  daughterInParentRestFrame(float parentLambdaPx, float parentLambdaPy, float parentLambdaPz, float parentLambdaMass,  float daughterProtonPx, float daughterProtonPy, float daughterProtonPz , float daughterProtonMass ){
      
        ROOT::Math::PxPyPzMVector parentFourVector{parentLambdaPx, parentLambdaPy, parentLambdaPz, parentLambdaMass};
        ROOT::Math::PxPyPzMVector daughterFourVector{daughterProtonPx, daughterProtonPy, daughterProtonPz, daughterProtonMass};

        ROOT::Math::Boost boostToParentRest{parentFourVector.BoostToCM()}; //returns inverse boost requiret make lambda momentum zero i.e rest fram of lambda
        const auto protonStar = boostToParentRest(daughterFourVector); //Applies the store inverse boost of the lambda rest fram to the proton
      
        return protonStar;
    } 

    float CosTheta (ROOT::Math::PxPyPzMVector proton1, ROOT::Math::PxPyPzMVector proton2){
      float numerator =  ( proton1.Px()*proton2.Px() + proton1.Py()*proton2.Py() + proton1.Pz()*proton2.Pz() );
      float denominator = (proton1.P()*proton2.P() );
      float cosTheta;
      if (denominator == 0.f){
        cosTheta = -2;
      }else{cosTheta = numerator/denominator;}
      return cosTheta;
    }

  template <typename Candidate1, typename Candidate2 >
  bool isKinematicallyCompatible(Candidate1 const& candidate1, Candidate2 const& candidate2)
  {
    const float deltaPt = std::abs(candidate1.pt() -candidate2.pt());
    const float deltaPhi = std::abs(std::remainder( candidate1.phi() - candidate2.phi(), 2.f * M_PI));
    const float deltaRapidity = std::abs(candidate1.rapidity() - candidate2.rapidity());
    
    return deltaPt < compatibilityDeltaPt &&
       deltaPhi < compatibilityDeltaPhi &&
       deltaRapidity < compatibilityDeltaRapidity;  
  }

  template <typename LambdaSlice, typename AntiLambdaSlice>
  void fillLambdaAntiLambdaSameEvent(LambdaSlice const& lambdas, AntiLambdaSlice const& antiLambdas)
  { 
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/EventCount"), 0 );
      
      //Lambda-AntiLambda pairs
      for (auto const& [lambda, antiLambda] : combinations(CombinationsFullIndexPolicy( lambdas, antiLambdas))) {
      
      // Boost the proton into the Lambda rest frame.
      const auto protonStar = daughterInParentRestFrame( lambda.px(), lambda.py(), lambda.pz(), lambda.mass(), lambda.protonPx(), lambda.protonPy(), lambda.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar = daughterInParentRestFrame( antiLambda.px(), antiLambda.py(), antiLambda.pz(), antiLambda.mass(), antiLambda.protonPx(), antiLambda.protonPy(), antiLambda.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar, antiProtonStar);
      const float deltaPhi = std::remainder(lambda.phi() - antiLambda.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/NumberOfPairs"), 0 );
      
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaEtaDelatPhi"), lambda.eta() - antiLambda.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaPtDelatEta"), lambda.pt() - antiLambda.pt(), lambda.eta() - antiLambda.eta() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaPtDelatPhi"), lambda.pt() - antiLambda.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaMassDeltaCent"), lambda.mass() - antiLambda.mass(), lambda.centrality() - antiLambda.centrality() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaMassDeltaPt"), lambda.mass() - antiLambda.mass(), lambda.pt() - antiLambda.pt() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/DeltaCentDeltaMult"), lambda.centrality() - antiLambda.centrality(), lambda.multiplicity() - antiLambda.multiplicity() );

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/Mass"), lambda.mass(), antiLambda.mass());

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/CosMult"), cosDeltaThetaStar, lambda.multiplicity());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/CosCent"), cosDeltaThetaStar, lambda.centrality());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hSameEvent/Cos"), cosDeltaThetaStar);
      }

    }

  template <typename LambdaSlice>
  void fillLambdaLambdaSameEvent(LambdaSlice const& lambdas)
  { 
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/EventCount"), 0 );
      
      //Lambda-Lamda pairs
      for (auto const& [lambda1, lambda2] : combinations(CombinationsStrictlyUpperIndexPolicy( lambdas, lambdas))) {
      
      // Boost the proton into the Lambda rest frame.
      const auto protonStar1 = daughterInParentRestFrame( lambda1.px(), lambda1.py(), lambda1.pz(), lambda1.mass(), lambda1.protonPx(), lambda1.protonPy(), lambda1.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto protonStar2 = daughterInParentRestFrame( lambda2.px(), lambda2.py(), lambda2.pz(), lambda2.mass(), lambda2.protonPx(), lambda2.protonPy(), lambda2.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar1, protonStar2);
      const float deltaPhi = std::remainder(lambda1.phi() - lambda2.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}

      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/NumberOfPairs"), 0 );
      
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaEtaDelatPhi"), lambda1.eta() - lambda2.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaPtDelatEta"), lambda1.pt() - lambda2.pt(), lambda1.eta() - lambda2.eta() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaPtDelatPhi"), lambda1.pt() - lambda2.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaMassDeltaCent"), lambda1.mass() - lambda2.mass(), lambda1.centrality() - lambda2.centrality() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaMassDeltaPt"), lambda1.mass() - lambda2.mass(), lambda1.pt() - lambda2.pt() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/DeltaCentDeltaMult"), lambda1.centrality() - lambda2.centrality(), lambda1.multiplicity() - lambda2.multiplicity() );

      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/Mass"), lambda1.mass(), lambda2.mass());

      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/CosMult"), cosDeltaThetaStar, lambda1.multiplicity());
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/CosCent"), cosDeltaThetaStar, lambda1.centrality());
      rSpinAnalysis.fill(HIST( "LambdaLambda/hSameEvent/Cos"), cosDeltaThetaStar);
      }

    }

  template <typename AntiLambdaSlice>
  void fillAntiLambdaAntiLambdaSameEvent( AntiLambdaSlice const& antiLambdas)
  { 
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/EventCount"), 0 );
       
      //AntiLambda-AntiLambda pairs
      for (auto const& [antilambda1, antilambda2] : combinations(CombinationsStrictlyUpperIndexPolicy( antiLambdas, antiLambdas))) {

      // Boost the proton into the Lambda rest frame.
      const auto antiProtonStar1 = daughterInParentRestFrame( antilambda1.px(), antilambda1.py(), antilambda1.pz(), antilambda1.mass(), antilambda1.protonPx(), antilambda1.protonPy(), antilambda1.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar2 = daughterInParentRestFrame( antilambda2.px(), antilambda2.py(), antilambda2.pz(), antilambda2.mass(), antilambda2.protonPx(), antilambda2.protonPy(), antilambda2.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( antiProtonStar1, antiProtonStar2);
      const float deltaPhi = std::remainder(antilambda1.phi() - antilambda2.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/NumberOfPairs"), 0 );

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaEtaDelatPhi"), antilambda1.eta() - antilambda2.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaPtDelatEta"), antilambda1.pt() - antilambda2.pt(), antilambda1.eta() - antilambda2.eta() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaPtDelatPhi"), antilambda1.pt() - antilambda2.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaMassDeltaCent"), antilambda1.mass() - antilambda2.mass(), antilambda1.centrality() - antilambda2.centrality() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaMassDeltaPt"), antilambda1.mass() - antilambda2.mass(), antilambda1.pt() - antilambda2.pt() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/DeltaCentDeltaMult"), antilambda1.centrality() - antilambda2.centrality(), antilambda1.multiplicity() - antilambda2.multiplicity() );

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/Mass"), antilambda1.mass(), antilambda2.mass());

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/CosMult"), cosDeltaThetaStar, antilambda1.multiplicity());
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/CosCent"), cosDeltaThetaStar, antilambda1.centrality());
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hSameEvent/Cos"), cosDeltaThetaStar);
      }
    }

  template <typename LambdaSlice1, typename AntiLambdaSlice2>
  void fillLambdaAntiLambdaMixedEvent(LambdaSlice1 const& lambdas, AntiLambdaSlice2 const& antiLambdas)
  { 
      bool hasAcceptedPair = false;
      //Lambda-AntiLambda pairs
      for (auto const& [lambda, antiLambda] : combinations(CombinationsFullIndexPolicy( lambdas, antiLambdas))) {
      
      if( !isKinematicallyCompatible(lambda, antiLambda) ){continue;}
      // Boost the proton into the Lambda rest frame.
      const auto protonStar = daughterInParentRestFrame( lambda.px(), lambda.py(), lambda.pz(), lambda.mass(), lambda.protonPx(), lambda.protonPy(), lambda.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar = daughterInParentRestFrame( antiLambda.px(), antiLambda.py(), antiLambda.pz(), antiLambda.mass(), antiLambda.protonPx(), antiLambda.protonPy(), antiLambda.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar, antiProtonStar);
      const float deltaPhi = std::remainder(lambda.phi() - antiLambda.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}
      hasAcceptedPair = true;

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/NumberOfPairs"), 0 );

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaEtaDelatPhi"), lambda.eta() - antiLambda.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaPtDelatEta"), lambda.pt() - antiLambda.pt(), lambda.eta() - antiLambda.eta() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaPtDelatPhi"), lambda.pt() - antiLambda.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaMassDeltaCent"), lambda.mass() - antiLambda.mass(), lambda.centrality() - antiLambda.centrality() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaMassDeltaPt"), lambda.mass() - antiLambda.mass(), lambda.pt() - antiLambda.pt() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaCentDeltaMult"), lambda.centrality() - antiLambda.centrality(), lambda.multiplicity() - antiLambda.multiplicity() );

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/Mass"), lambda.mass(), antiLambda.mass());

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/CosMult"), cosDeltaThetaStar, lambda.multiplicity(), antiLambda.multiplicity());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/CosCent"), cosDeltaThetaStar, lambda.centrality(), antiLambda.centrality());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/Cos"), cosDeltaThetaStar);
      }

      if (hasAcceptedPair) {
        rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/EventCount"), 0 );
      }

    }


  template < typename AntiLambdaSlice1, typename LambdaSlice2>
  void fillAntiLambdaLambdaMixedEvent(AntiLambdaSlice1 const& antiLambdas, LambdaSlice2 const& lambdas)
  { 
      bool hasAcceptedPair = false;
      //Lambda-AntiLambda pairs
      for (auto const& [lambda, antiLambda] : combinations(CombinationsFullIndexPolicy( lambdas, antiLambdas))) {
      
      if( !isKinematicallyCompatible(lambda, antiLambda) ){continue;}
      // Boost the proton into the Lambda rest frame.
      const auto protonStar = daughterInParentRestFrame( lambda.px(), lambda.py(), lambda.pz(), lambda.mass(), lambda.protonPx(), lambda.protonPy(), lambda.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar = daughterInParentRestFrame( antiLambda.px(), antiLambda.py(), antiLambda.pz(), antiLambda.mass(), antiLambda.protonPx(), antiLambda.protonPy(), antiLambda.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar, antiProtonStar);
      const float deltaPhi = std::remainder(lambda.phi() - antiLambda.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}
      hasAcceptedPair = true;

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/NumberOfPairs"), 0 );
      
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaEtaDelatPhi"), lambda.eta() - antiLambda.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaPtDelatEta"), lambda.pt() - antiLambda.pt(), lambda.eta() - antiLambda.eta() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaPtDelatPhi"), lambda.pt() - antiLambda.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaMassDeltaCent"), lambda.mass() - antiLambda.mass(), lambda.centrality() - antiLambda.centrality() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaMassDeltaPt"), lambda.mass() - antiLambda.mass(), lambda.pt() - antiLambda.pt() );
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/DeltaCentDeltaMult"), lambda.centrality() - antiLambda.centrality(), lambda.multiplicity() - antiLambda.multiplicity() );

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/Mass"), lambda.mass(), antiLambda.mass());

      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/CosMult"), cosDeltaThetaStar, lambda.multiplicity(), antiLambda.multiplicity());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/CosCent"), cosDeltaThetaStar, lambda.centrality(), antiLambda.centrality());
      rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/Cos"), cosDeltaThetaStar);
      }

      if (hasAcceptedPair) {
        rSpinAnalysis.fill(HIST( "LambdaAntiLambda/hMixedEvent/EventCount"), 0 );
      }
    }


  template <typename LambdaSlice1, typename LambdaSlice2>
  void fillLambdaLambdaMixedEvent(LambdaSlice1 const& lambdas1, LambdaSlice2 const& lambdas2)
  { 
      bool hasAcceptedPair = false;
      //Lambda-AntiLambda pairs
      for (auto const& [lambda1, lambda2] : combinations(CombinationsFullIndexPolicy( lambdas1, lambdas2))) {
      
      if( !isKinematicallyCompatible(lambda1, lambda2) ){continue;}
      // Boost the proton into the Lambda rest frame.
      const auto protonStar = daughterInParentRestFrame( lambda1.px(), lambda1.py(), lambda1.pz(), lambda1.mass(), lambda1.protonPx(), lambda1.protonPy(), lambda1.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar = daughterInParentRestFrame( lambda2.px(), lambda2.py(), lambda2.pz(), lambda2.mass(), lambda2.protonPx(), lambda2.protonPy(), lambda2.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar, antiProtonStar);
      const float deltaPhi = std::remainder(lambda1.phi() - lambda2.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}
      hasAcceptedPair = true;

      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/NumberOfPairs"), 0 );

      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaEtaDelatPhi"), lambda1.eta() - lambda2.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaPtDelatEta"), lambda1.pt() - lambda2.pt(), lambda1.eta() - lambda2.eta() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaPtDelatPhi"), lambda1.pt() - lambda2.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaMassDeltaCent"), lambda1.mass() - lambda2.mass(), lambda1.centrality() - lambda2.centrality() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaMassDeltaPt"), lambda1.mass() - lambda2.mass(), lambda1.pt() - lambda2.pt() );
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/DeltaCentDeltaMult"), lambda1.centrality() - lambda2.centrality(), lambda1.multiplicity() - lambda2.multiplicity() );

      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/Mass"), lambda1.mass(), lambda2.mass());

      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/CosMult"), cosDeltaThetaStar, lambda1.multiplicity(), lambda2.multiplicity());
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/CosCent"), cosDeltaThetaStar, lambda1.centrality(), lambda2.centrality());
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/Cos"), cosDeltaThetaStar);
      }

      if (hasAcceptedPair) {
      rSpinAnalysis.fill(HIST( "LambdaLambda/hMixedEvent/EventCount"), 0 );
      }

    }

  template <typename AntiLambdaSlice1, typename AntiLambdaSlice2>
  void fillAntiLambdaAntiLambdaMixedEvent(AntiLambdaSlice1 const& antiLambdas1, AntiLambdaSlice2 const& antiLambdas2)
  { 
      bool hasAcceptedPair = false;
      //Lambda-AntiLambda pairs
      for (auto const& [antiLambda1, antiLambda2] : combinations(CombinationsFullIndexPolicy( antiLambdas1, antiLambdas2))) {
      
      if( !isKinematicallyCompatible(antiLambda1, antiLambda2) ){continue;}
      // Boost the proton into the Lambda rest frame.
      const auto protonStar = daughterInParentRestFrame( antiLambda1.px(), antiLambda1.py(), antiLambda1.pz(), antiLambda1.mass(), antiLambda1.protonPx(), antiLambda1.protonPy(), antiLambda1.protonPz(), o2::constants::physics::MassProton);
      // Boost the antiproton into the anti-Lambda rest frame.
      const auto antiProtonStar = daughterInParentRestFrame( antiLambda2.px(), antiLambda2.py(), antiLambda2.pz(), antiLambda2.mass(), antiLambda2.protonPx(), antiLambda2.protonPy(), antiLambda2.protonPz(), o2::constants::physics::MassProton);

      const float cosDeltaThetaStar = CosTheta( protonStar, antiProtonStar);
      const float deltaPhi = std::remainder(antiLambda1.phi() - antiLambda2.phi(), 2.f * M_PI);
      // Invalid result returned by cosTheta().
      if (cosDeltaThetaStar < -1.f) {continue;}
      hasAcceptedPair = true;

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/NumberOfPairs"), 0 );

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaEtaDelatPhi"), antiLambda1.eta() - antiLambda2.eta(), deltaPhi );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaPtDelatEta"), antiLambda1.pt() - antiLambda2.pt(), antiLambda1.eta() - antiLambda2.eta() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaPtDelatPhi"), antiLambda1.pt() - antiLambda2.pt(), deltaPhi );


      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaMassDeltaCent"), antiLambda1.mass() - antiLambda2.mass(), antiLambda1.centrality() - antiLambda2.centrality() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaMassDeltaPt"), antiLambda1.mass() - antiLambda2.mass(), antiLambda1.pt() - antiLambda2.pt() );
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/DeltaCentDeltaMult"), antiLambda1.centrality() - antiLambda2.centrality(), antiLambda1.multiplicity() - antiLambda2.multiplicity() );

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/Mass"), antiLambda1.mass(), antiLambda2.mass());

      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/CosMult"), cosDeltaThetaStar, antiLambda1.multiplicity(), antiLambda2.multiplicity());
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/CosCent"), cosDeltaThetaStar, antiLambda1.centrality(), antiLambda2.centrality());
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/Cos"), cosDeltaThetaStar);
      }

      if (hasAcceptedPair) {
      rSpinAnalysis.fill(HIST( "AntiLambdaAntiLambda/hMixedEvent/EventCount"), 0 );
      }

    }
  // ================================================================
  // Mass Filter
  // ================================================================

  Configurable<float> minLambdaMass{"minLambdaMass", 1.11f, "Minimum Lambda mass"};
  Configurable<float> maxLambdaMass{"maxLambdaMass", 1.12f, "Maximum Lambda mass"};


  Filter lambdaMassFilter = aod::lambdahyperon::mass >= minLambdaMass &&
                            aod::lambdahyperon::mass <= maxLambdaMass;

  using FilteredLambdas = soa::Filtered<aod::Lambdas>;
  using FilteredAntiLambdas = soa::Filtered<aod::AntiLambdas>;

  Preslice<FilteredLambdas> lambdasPerCollision = aod::lambdahyperon::collisionId;
  Preslice<FilteredAntiLambdas> antiLambdasPerCollision = aod::lambdahyperon::collisionId;

    void process( aod::Collisions const& collisions, FilteredLambdas const& lambdas, FilteredAntiLambdas const& antiLambdas){
      
      for (auto const& collision : collisions) {

        const auto lambdasThisCollision = lambdas.sliceBy( lambdasPerCollision, collision.globalIndex());

        const auto antiLambdasThisCollision = antiLambdas.sliceBy( antiLambdasPerCollision, collision.globalIndex());

        // Lambda–anti-Lambda requires at least one of each.
        if (lambdasThisCollision.size() > 0 && antiLambdasThisCollision.size() > 0) {
          fillLambdaAntiLambdaSameEvent(lambdasThisCollision, antiLambdasThisCollision);}

        // Lambda–Lambda requires at least two Lambdas.
        if (lambdasThisCollision.size() >= 2) {fillLambdaLambdaSameEvent(lambdasThisCollision);}

        // Anti-Lambda–anti-Lambda requires at least two anti-Lambdas.
        if (antiLambdasThisCollision.size() >= 2) {fillAntiLambdaAntiLambdaSameEvent(antiLambdasThisCollision);}
      
      }

      // ============================================================
      // Mixed events
      // ============================================================

      for ( auto const& [collision1, collision2] : combinations( CombinationsStrictlyUpperIndexPolicy(collisions, collisions) )  ) {

      const auto lambdas1 = lambdas.sliceBy(lambdasPerCollision, collision1.globalIndex());
      const auto lambdas2 = lambdas.sliceBy(lambdasPerCollision, collision2.globalIndex());

      const auto antiLambdas1 = antiLambdas.sliceBy(antiLambdasPerCollision, collision1.globalIndex());
      const auto antiLambdas2 = antiLambdas.sliceBy(antiLambdasPerCollision, collision2.globalIndex());

      // Lambda–anti-Lambda requires at least one of each.
        if ( lambdas1.size() > 0 && antiLambdas2.size() > 0 ) {
          fillLambdaAntiLambdaMixedEvent(lambdas1, antiLambdas2);
        }
        if ( lambdas2.size() > 0 && antiLambdas1.size() > 0 ){
          fillAntiLambdaLambdaMixedEvent(antiLambdas1, lambdas2);
        }
      // Lambda–Lambda requires at least two Lambdas.
        if (lambdas1.size() > 0 && lambdas2.size() > 0 ) {
          fillLambdaLambdaMixedEvent(lambdas1, lambdas2);
        }
      // Anti-Lambda–anti-Lambda requires at least two anti-Lambdas.
        if (antiLambdas1.size() > 0 && antiLambdas2.size() > 0 ) {
          fillAntiLambdaAntiLambdaMixedEvent(antiLambdas1, antiLambdas2);
        }
      }
  }

};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  const bool isMC = cfgc.options().get<bool>("isMC");

  WorkflowSpec workflow;

  if (isMC) {
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaMcRecoTableProducer>(cfgc));
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaEfficiencyPlots>(cfgc));
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaSelectionCutFlow>(cfgc));
  } else {
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaProducer>(cfgc));
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaSelector>(cfgc));
    workflow.push_back( adaptAnalysisTask<LambdaAntiLambdaPairAnalysis>(cfgc));
  }

  return workflow;
}

/*WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LambdaAntiLambdaProducer>(cfgc),
    adaptAnalysisTask<LambdaAntiLambdaSelector>(cfgc),
    adaptAnalysisTask<LambdaAntiLambdaPairAnalysis>(cfgc)};
    //adaptAnalysisTask<LambdaAntiLambdaMcRecoTableProducer>(cfgc),
    //adaptAnalysisTask<LambdaAntiLambdaEfficiencyPlots>(cfgc),
    //adaptAnalysisTask<LambdaAntiLambdaSelectionCutFlow>(cfgc)
}*/
