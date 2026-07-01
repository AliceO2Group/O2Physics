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

/// \file   lambda1405analysis.cxx
/// \brief Analysis task for lambda1405 via sigma kink decay
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFLambda1405Table.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>

#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullPr>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullWCentQVecs = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::QvectorFT0Cs>;
using CollisionsFullMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
using CollisionsFullMcWCent = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;

enum L1405DecayType { kSigmaMinusPiToPiPiNeutron = 0,
                      kSigmaPlusPiToPiPiNeutron,
                      kSigmaPlusPiToPiPiProton,
                      kNL1405Decays };

struct lambda1405candidate {
  // Columns for Lambda(1405) candidate
  float massL1405 = -1;                                     // Invariant mass of the Lambda(1405) candidate
  float massXi1530 = -1;                                    // Invariant mass of the Xi(1530) candidate
  float px = -1;                                            // Px of the Lambda(1405) candidate
  float py = -1;                                            // Py of the Lambda(1405) candidate
  float pz = -1;                                            // Pz of the Lambda(1405) candidate
  float pt() const { return std::sqrt(px * px + py * py); } // pT of the Lambda(1405) candidate
  float phi = -1;                                           // Phi of the Lambda(1405) candidate

  bool isSigmaPlus = false;   // True if compatible with Sigma+
  bool isSigmaMinus = false;  // True if compatible with Sigma-
  float sigmaMinusMass = -1;  // Invariant mass of the Sigma- candidate
  float sigmaPlusMass = -1;   // Invariant mass of the Sigma+ candidate
  float xiMinusMass = -1;     // Invariant mass of the Xi- candidate
  int sigmaSign = 0;          // Sign of the Sigma candidate: 1 for matter, -1 for antimatter
  float sigmaPt = -1;         // pT of the Sigma daughter
  float sigmaAlphaAP = -1;    // Alpha of the Sigma
  float sigmaQtAP = -1;       // qT of the Sigma
  float sigmaRadius = -1;     // Radius of the Sigma decay vertex
  float dcaSigmaToPv = -1;    // DCA of the Sigma candidate to the primary vertex
  float kinkPt = -1;          // pT of the kink daughter
  float kinkPiNSigTpc = -1;   // Number of sigmas for the pion candidate from Sigma kink in Tpc
  float kinkPiNSigTof = -1;   // Number of sigmas for the pion candidate from Sigma kink in Tof
  float kinkPrNSigTpc = -1;   // Number of sigmas for the proton candidate from Sigma kink in Tpc
  float kinkPrNSigTof = -1;   // Number of sigmas for the proton candidate from Sigma kink in Tof
  float kinkDcaDauToPv = -1;  // DCA of the kink daughter to the primary vertex

  float bachPiPt = -1;      // pT of the pion daughter
  float bachPiNSigTpc = -1; // Number of sigmas for the pion candidate
  float bachPiNSigTof = -1; // Number of sigmas for the pion candidate using Tof
  int kinkDauId = 0;        // Id of the pion from Sigma decay in MC
  int sigmaId = 0;          // Id of the Sigma candidate in MC
  int bachPiId = 0;         // Id of the pion candidate in MC

  float cent = -1;      // Centrality of the collision
  float pvContrib = -1; // Number of contributors to the primary vertex

  float scalarProd = -1; // Scalar product for flow analysis
};

struct lambda1405analysis {
  int lambda1405PdgCode = 102132;                     // PDG code for Lambda(1405)
  int sigmaMinusPdgCode = 3112;                       // PDG code for Lambda(1405)
  int sigmaPlusPdgCode = 3222;                        // PDG code for Lambda(1405)
  Produces<aod::Lambda1405Cands> outputDataTable;     // Output table for Lambda(1405) candidates
  Produces<aod::Lambda1405Flow> outputDataFlowTable;  // Output table for Lambda(1405) flow analysis
  Produces<aod::Lambda1405CandsMC> outputDataTableMC; // Output table for Lambda(1405) candidates in MC
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda1405{"lambda1405", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaMinus{"sigmaMinus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaPlus{"sigmaPlus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSelections{"selections", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutEtaDaughter{"cutEtaDaughter", 0.8f, "Eta cut for daughter tracks"};
  Configurable<float> cutDCAtoPvSigma{"cutDCAtoPvSigma", 0.1f, "Max DCA to primary vertex for Sigma candidates (cm)"};
  Configurable<float> cutDCAtoPvPiFromSigma{"cutDCAtoPvPiFromSigma", 2., "Min DCA to primary vertex for pion from Sigma candidates (cm)"};

  Configurable<float> cutMinPtL1405{"cutMinPtL1405", 2.0f, "Minimum pT cut for Lambda(1405) candidates (GeV/c)"};
  Configurable<float> cutUpperMass{"cutUpperMass", 1.6f, "Upper mass cut for Lambda(1405) candidates (GeV/c^2)"};
  Configurable<float> cutSigmaRadius{"cutSigmaRadius", 20.f, "Minimum radius for Sigma candidates (cm)"};
  Configurable<float> cutSigmaMass{"cutSigmaMass", 0.1, "Sigma mass window (MeV/c^2)"};
  Configurable<float> cutNITSClusKink{"cutNITSClusKink", 3, "Minimum number of ITS clusters for pion candidate"};
  Configurable<float> cutNTpcClusPi{"cutNTpcClusPi", 90, "Minimum number of Tpc clusters for pion candidate"};
  Configurable<float> cutNSigTpc{"cutNSigTpc", 3, "NSigTpcPion"};
  Configurable<float> cutNSigTof{"cutNSigTof", 3, "NSigTofPion"};
  Configurable<std::string> cutSigmaQtAPMin{"cutSigmaQtAPMin", "0.17", "Functional form of minimum qT(alpha) selection"};
  Configurable<std::string> cutSigmaQtAPMax{"cutSigmaQtAPMax", "0.2", "Functional form of minimum qT(alpha) selection"};
  Configurable<std::string> cutMinBachPiPtVsL1405Pt{"cutMinBachPiPtVsL1405Pt", "0", "Functional form of minimum qT(alpha) selection"};
  Configurable<std::string> cutMaxBachPiPtVsL1405Pt{"cutMaxBachPiPtVsL1405Pt", "10", "Functional form of minimum qT(alpha) selection"};
  Configurable<std::string> cutMinSigmaPtVsL1405Pt{"cutMinSigmaPtVsL1405Pt", "0", "Functional form of minimum qT(alpha) selection"};
  Configurable<std::string> cutMaxSigmaPtVsL1405Pt{"cutMaxSigmaPtVsL1405Pt", "10", "Functional form of minimum qT(alpha) selection"};

  // Configurables for flow analysis
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP/EP computation (not applied in resolution process)"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP/EP computation (not applied in resolution process)"};
  Configurable<bool> fillFlowTree{"fillFlowTree", false, "If true, fill the output tree with Lambda(1405) candidates"};

  // Downsample for table filling
  Configurable<float> downSampleFactor{"downSampleFactor", 1., "Fraction of candidates to keep in TTree"};
  Configurable<float> ptDownSampleMax{"ptDownSampleMax", 10., "Maximum pt for the application of the downsampling factor"};

  Configurable<bool> fillOutputTree{"fillOutputTree", true, "If true, fill the output tree with Lambda(1405) candidates"};
  Configurable<bool> doLikeSignBkg{"doLikeSignBkg", false, "Use like-sign background"};
  Configurable<bool> useTof{"useTof", false, "Use Tof for PId for pion candidates"};

  TF1 funcMinQtAlphaAP, funcMaxQtAlphaAP, funcMinBachPiPtVsL1405Pt, funcMaxBachPiPtVsL1405Pt, funcMinSigmaPtVsL1405Pt, funcMaxSigmaPtVsL1405Pt;

  using CollisionsCentSel = soa::Filtered<CollisionsFullWCentQVecs>;
  using McRecoCollisionsCentSel = soa::Filtered<CollisionsFullMcWCent>;

  Filter filterCentrality = aod::cent::centFT0C >= centralityMin && aod::cent::centFT0C <= centralityMax;

  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // Configurable axes
  ConfigurableAxis axisPvContrib{"axisPvContrib", {5000, 0., 5000.}, ""};
  ConfigurableAxis axisCent{"axisCent", {100, 0., 100.}, ""};
  ConfigurableAxis axisPtL1405{"axisPtL1405", {100, 0., 10.}, ""};
  ConfigurableAxis axisPCompsL1405{"axisPCompsL1405", {100, -10., 10.}, ""};
  ConfigurableAxis axisPtKinkDaug{"axisPtKinkDaug", {100, 0., 10.}, ""};
  ConfigurableAxis axisPtResolution{"axisPtResolution", {100, -0.5, 0.5}, ""};
  ConfigurableAxis axisMassL1405{"axisMassL1405", {100, 1.3, 1.8}, ""};
  ConfigurableAxis axisXi1530Mass{"axisXi1530Mass", {100, 1.4, 1.6}, ""};
  ConfigurableAxis axisMassResolution{"axisMassResolution", {100, -0.1, 0.1}, ""};
  ConfigurableAxis axisNSig{"axisNSig", {100, -5., 5.}, ""};
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {100, 1.1, 1.4}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {100, 1.2, 1.5}, ""};
  ConfigurableAxis axisVertexZ{"axisVertexZ", {100, -15., 15.}, ""};
  ConfigurableAxis axisQtAP{"axisQtAP", {100, 0., 0.3}, ""};
  ConfigurableAxis axisAlphaAP{"axisAlphaAP", {200, -1., 1.}, ""};
  ConfigurableAxis axisSigmaRadius{"axisSigmaRadius", {100, 0., 100.}, ""};
  ConfigurableAxis axisDcaSigmaToPv{"axisDcaSigmaToPv", {120, 0., 0.05}, ""};
  ConfigurableAxis axisDcaKinkToPv{"axisDcaKinkToPv", {200, 0., 20.}, ""};
  ConfigurableAxis axisScalarProd{"axisScalarProd", {200, -4., 4.}, ""};

  Preslice<aod::KinkCands> mKinkPerCol = aod::track::collisionId;
  Preslice<aod::TracksIU> mPerColTracks = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{axisPtL1405, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec pCompAxisL1405{axisPCompsL1405, "#it{p}_{comp} (GeV/#it{c})"};
    const AxisSpec ptKinkDaugAxis{axisPtKinkDaug, "#it{p}_{T}^{#pi} (GeV/#it{c})"};
    const AxisSpec ptResolutionAxis{axisPtResolution, "#it{p}_{T}^{rec} - #it{p}_{T}^{gen} (GeV/#it{c})"};
    const AxisSpec lambda1405MassAxis{axisMassL1405, "m(#Sigma#pi) (GeV/#it{c}^{2})"};
    const AxisSpec xi1530MassAxis{axisXi1530Mass, "m(#Xi#pi) (GeV/#it{c}^{2})"};
    const AxisSpec massResolutionAxis{axisMassResolution, "m_{rec} - m_{gen} (GeV/#it{c}^{2})"};
    const AxisSpec nSigmaAxis{axisNSig, "n#sigma_{#pi}"};
    const AxisSpec sigmaMassAxis{axisSigmaMass, "m(#Sigma) (GeV/#it{c}^{2})"};
    const AxisSpec xiMassAxis{axisXiMass, "m(#Xi) (GeV/#it{c}^{2})"};
    const AxisSpec vertexZAxis{axisVertexZ, "vtx_{Z} (cm)"};
    const AxisSpec qtAxis{axisQtAP, "q_{T, AP}"};
    const AxisSpec alphaAxis{axisAlphaAP, "#alpha_{AP}"};
    const AxisSpec sigmaRadiusAxis{axisSigmaRadius, "Dec. radius #Sigma (cm)"};
    const AxisSpec centAxis{axisCent, "Centrality (%)"};
    const AxisSpec pvContribAxis{axisPvContrib, "PV Contributors"};
    const AxisSpec dcaSigmaToPvBinsAxis{axisDcaSigmaToPv, "#Sigma DCA to PV (cm)"};
    const AxisSpec dcaKinkToPvBinsAxis{axisDcaKinkToPv, "Kink daug. DCA to PV (cm)"};
    const AxisSpec scalarProdAxis{axisScalarProd, "SP"};

    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // Sigma- candidate properties
    rSigmaMinus.add("hMassXiMinusSigmaMinus", "hMassXiMinusSigmaMinus", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2PtMassSigmaMinusBeforeCuts", "h2PtMassSigmaMinusBeforeCuts", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2PtPiKinkNSigBeforeCutsSigmaMinus", "h2PtPiKinkNSigBeforeCutsSigmaMinus", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
    rSigmaMinus.add("h2dPtMassSigmaMinus", "h2dPtMassSigmaMinus", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2dPvContribMassSigmaMinus", "h2dPvContribMassSigmaMinus", {HistType::kTH2F, {pvContribAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2dPvContribPtSigmaMinus", "h2dPvContribPtSigmaMinus", {HistType::kTH2F, {pvContribAxis, ptAxis}});
    rSigmaMinus.add("h2KinkPiPtNSigTofSigmaMinus", "h2KinkPiPtNSigTofSigmaMinus", {HistType::kTH2F, {ptKinkDaugAxis, nSigmaAxis}});
    rSigmaMinus.add("hSigmaMinusArmPod", "hSigmaMinusArmPod", {HistType::kTH2D, {alphaAxis, qtAxis}});
    rSigmaMinus.add("hSigmaMinusRadius", "hSigmaMinusRadius", {HistType::kTH2F, {sigmaRadiusAxis, ptAxis}});
    rSigmaMinus.add("hSigmaMinusDcaToPv", "hSigmaMinusDcaToPv", {HistType::kTH2F, {dcaSigmaToPvBinsAxis, ptAxis}});
    rSigmaMinus.add("hSigmaMinusKinkTpcNSigPi", "hSigmaMinusKinkTpcNSigPi", {HistType::kTH2F, {nSigmaAxis, ptKinkDaugAxis}});
    rSigmaMinus.add("hSigmaMinusKinkTofNSigPi", "hSigmaMinusKinkTofNSigPi", {HistType::kTH2F, {nSigmaAxis, ptKinkDaugAxis}});
    rSigmaMinus.add("hSigmaMinusDcaKinkDauToPv", "hSigmaMinusDcaKinkDauToPv", {HistType::kTH2F, {dcaKinkToPvBinsAxis, ptAxis}});

    // Sigma+ candidate properties
    rSigmaPlus.add("h2PtMassSigmaPlusBeforeCuts", "h2PtMassSigmaPlusBeforeCuts", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaPlus.add("h2PtPrKinkNSigBeforeCuts", "h2PtPrKinkNSigBeforeCuts", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
    rSigmaPlus.add("h2dPtMassSigmaPlus", "h2dPtMassSigmaPlus", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaPlus.add("h2dPvContribMassSigmaPlus", "h2dPvContribMassSigmaPlus", {HistType::kTH2F, {pvContribAxis, sigmaMassAxis}});
    rSigmaPlus.add("h2dPvContribPtSigmaPlus", "h2dPvContribPtSigmaPlus", {HistType::kTH2F, {pvContribAxis, ptAxis}});
    rSigmaPlus.add("h2KinkPrPtNSigTofSigmaPlus", "h2KinkPrPtNSigTofSigmaPlus", {HistType::kTH2F, {ptKinkDaugAxis, nSigmaAxis}});
    rSigmaPlus.add("hSigmaPlusArmPod", "hSigmaPlusArmPod", {HistType::kTH2F, {alphaAxis, qtAxis}});
    rSigmaPlus.add("hSigmaPlusRadius", "hSigmaPlusRadius", {HistType::kTH2F, {sigmaRadiusAxis, ptAxis}});
    rSigmaPlus.add("hSigmaPlusDcaToPv", "hSigmaPlusDcaToPv", {HistType::kTH2F, {dcaSigmaToPvBinsAxis, ptAxis}});
    rSigmaPlus.add("hSigmaPlusKinkTpcNSigPr", "hSigmaPlusKinkTpcNSigPr", {HistType::kTH2F, {nSigmaAxis, ptKinkDaugAxis}});
    rSigmaPlus.add("hSigmaPlusKinkTofNSigPr", "hSigmaPlusKinkTofNSigPr", {HistType::kTH2F, {nSigmaAxis, ptKinkDaugAxis}});
    rSigmaPlus.add("hSigmaPlusDcaKinkDauToPv", "hSigmaPlusDcaKinkDauToPv", {HistType::kTH2F, {dcaKinkToPvBinsAxis, ptKinkDaugAxis}});
    rSigmaPlus.add("hMassXiMinusSigmaPlus", "hMassXiMinusSigmaPlus", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});

    // Mass QA
    rLambda1405.add("h2SigmaMinusMassVsLambdaMass", "h2SigmaMinusMassVsLambdaMass", {HistType::kTH2F, {lambda1405MassAxis, sigmaMassAxis}});
    rLambda1405.add("h2SigmaPlusMassVsLambdaMass", "h2SigmaPlusMassVsLambdaMass", {HistType::kTH2F, {lambda1405MassAxis, sigmaMassAxis}});
    rLambda1405.add("hMassL1405", "hMassL1405", {HistType::kTH2F, {lambda1405MassAxis, ptAxis}});
    rLambda1405.add("hMassXi1530", "hMassXi1530", {HistType::kTH2F, {xi1530MassAxis, ptAxis}});
    // Kinematic distributions
    rLambda1405.add("hPx", "hPx;#it{p}_x;Counts", {HistType::kTH1D, {pCompAxisL1405}});
    rLambda1405.add("hPy", "hPy;#it{p}_y;Counts", {HistType::kTH1D, {pCompAxisL1405}});
    rLambda1405.add("hPz", "hPz;#it{p}_z;Counts", {HistType::kTH1D, {pCompAxisL1405}});
    rLambda1405.add("hPt", "hPt", {HistType::kTH1D, {ptAxis}});
    rLambda1405.add("hPhi", "hPhi", {HistType::kTH1D, {{128, -o2::constants::math::PI, o2::constants::math::PI}}});
    // Pion daughter properties
    rLambda1405.add("h2BachPiPtNSigTof", "h2BachPiPtNSigTof", {HistType::kTH2F, {ptKinkDaugAxis, nSigmaAxis}});
    rLambda1405.add("h2BachPiPtNSigTpc", "h2BachPiPtNSigTpc", {HistType::kTH2F, {ptKinkDaugAxis, nSigmaAxis}});
    // Sparse histogram
    std::vector<AxisSpec> axes = {lambda1405MassAxis, ptAxis, sigmaMassAxis, alphaAxis, qtAxis, dcaSigmaToPvBinsAxis, dcaKinkToPvBinsAxis};
    if (doprocessDataWCentQVecs) {
      axes.push_back(centAxis);
      axes.push_back(scalarProdAxis);
    }
    if (doprocessMc || doprocessMcWCentSel) {
      axes.push_back(pvContribAxis);
    }
    rLambda1405.add("hSparseL1405", "THn for SP", HistType::kTHnSparseF, axes);

    // Selection QA
    rSelections.add("hSelectionsL1405", "hSelectionsL1405", {HistType::kTH1D, {{6, -0.f, 5.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(2, "Passed Sigma sel");
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(3, "Passed Bach PID sel");
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(4, "Upper mass sel");
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(5, "#it{p}_{T} Correl");
    rSelections.get<TH1>(HIST("hSelectionsL1405"))->GetXaxis()->SetBinLabel(6, "Accepted");
    rSelections.add("hSelectionsSigmaPlus", "hSelectionsSigmaPlus", {HistType::kTH1D, {{7, -0.f, 6.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(2, "Passed kink sel");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(3, "Passed mass sel");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(4, "Passed cutDCAtoPvSigma");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(5, "Passed cutDCAtoPvPiFromSigma");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(6, "Passed cutSigmaRadius");
    rSelections.get<TH1>(HIST("hSelectionsSigmaPlus"))->GetXaxis()->SetBinLabel(7, "Passed cutSigmaQtAP");
    rSelections.add("hSelectionsSigmaMinus", "hSelectionsSigmaMinus", {HistType::kTH1D, {{7, -0.f, 6.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(2, "Passed kink sel");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(3, "Passed mass sel");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(4, "Passed cutDCAtoPvSigma");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(5, "Passed cutDCAtoPvPiFromSigma");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(6, "Passed cutSigmaRadius");
    rSelections.get<TH1>(HIST("hSelectionsSigmaMinus"))->GetXaxis()->SetBinLabel(7, "Passed cutSigmaQtAP");
    rSelections.add("hSelectionsBachPi", "hSelectionsBachPi", {HistType::kTH1D, {{6, -0.f, 5.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(2, "Sign sel");
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(3, "Nsigma Tpc sel");
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(4, "Tpc clusters sel");
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(5, "#eta sel");
    rSelections.get<TH1>(HIST("hSelectionsBachPi"))->GetXaxis()->SetBinLabel(6, "PID sel");
    rSelections.add("hSelectionsKinkPi", "hSelectionsKinkPi", {HistType::kTH1D, {{7, -0.f, 6.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(2, "Tpc N#sigma PID sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(3, "Tpc N clusters sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(4, "#eta sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(5, "ITS clusters sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(6, "has Tof sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPi"))->GetXaxis()->SetBinLabel(7, "N#sigma Tof sel");
    rSelections.add("hSelectionsKinkPr", "hSelectionsKinkPr", {HistType::kTH1D, {{7, -0.f, 6.5f}}});
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(1, "All");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(2, "Tpc N#sigma PID sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(3, "Tpc N clusters sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(4, "#eta sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(5, "ITS clusters sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(6, "has Tof sel");
    rSelections.get<TH1>(HIST("hSelectionsKinkPr"))->GetXaxis()->SetBinLabel(7, "N#sigma Tof sel");

    if (doprocessDataWCentQVecs) {
      rLambda1405.add("hScalarProd", "hScalarProd", {HistType::kTH2F, {scalarProdAxis, ptAxis}});
      std::vector<AxisSpec> axesFlow = {lambda1405MassAxis, ptAxis, centAxis, scalarProdAxis};
      rLambda1405.add("hSparseFlowL1405", "THn for SP", HistType::kTHnSparseF, axesFlow);
    }

    if (doprocessMc || doprocessMcWCentSel) {
      // Add MC histograms
      rSigmaPlus.add("h2GenSigmaPlusPvContribPt", "h2GenSigmaPlusPvContribPt", {HistType::kTH2F, {pvContribAxis, ptAxis}});
      rSigmaMinus.add("h2GenSigmaMinusPvContribPt", "h2GenSigmaMinusPvContribPt", {HistType::kTH2F, {pvContribAxis, ptAxis}});

      rLambda1405.add("hRecoL1405", "hRecoL1405;;Counts", {HistType::kTH2F, {{4, -0.5, 3.5}, ptAxis}});
      rLambda1405.get<TH2>(HIST("hRecoL1405"))->GetXaxis()->SetBinLabel(1, "Reco #Lambda(1405)");
      rLambda1405.get<TH2>(HIST("hRecoL1405"))->GetXaxis()->SetBinLabel(2, "Has bach #pi");
      rLambda1405.get<TH2>(HIST("hRecoL1405"))->GetXaxis()->SetBinLabel(3, "Has mothers");
      rLambda1405.get<TH2>(HIST("hRecoL1405"))->GetXaxis()->SetBinLabel(4, "Matched");
      rLambda1405.add("hGenL1405", "hGenL1405;;Counts", {HistType::kTH2F, {{3, -0.5, 2.5}, ptAxis}});
      rLambda1405.get<TH2>(HIST("hGenL1405"))->GetXaxis()->SetBinLabel(1, "#Lambda(1405) #rightarrow #Sigma^{-} #pi^{+} #rightarrow n #pi^{-} #pi^{+}");
      rLambda1405.get<TH2>(HIST("hGenL1405"))->GetXaxis()->SetBinLabel(2, "#Lambda(1405) #rightarrow #Sigma^{+} #pi^{-} #rightarrow n #pi^{+} #pi^{-}");
      rLambda1405.get<TH2>(HIST("hGenL1405"))->GetXaxis()->SetBinLabel(3, "#Lambda(1405) #rightarrow #Sigma^{+} #pi^{-} #rightarrow p #pi^{0} #pi^{-}");
      rLambda1405.add("h2GenL1405PvContribPt", "h2GenL1405PvContribPt", {HistType::kTH2F, {pvContribAxis, ptAxis}});
      rLambda1405.add("h2RecL1405PvContribPt", "h2RecL1405PvContribPt", {HistType::kTH2F, {pvContribAxis, ptAxis}});
      rLambda1405.add("h2GenSigmaMinusArmPod", "h2GenSigmaMinusArmPod", {HistType::kTH2F, {alphaAxis, qtAxis}});
      rLambda1405.add("h2GenSigmaPlusArmPod", "h2GenSigmaPlusArmPod", {HistType::kTH2F, {alphaAxis, qtAxis}});
      rLambda1405.add("h2GenPtVsBachPtSigmaMinus", "h2GenPtVsBachPtSigmaMinus;#Lambda(1405) #it{p}_{T} (GeV/c); Bach #pi #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenPtVsBachPtSigmaPlusToPi", "h2GenPtVsBachPtSigmaPlusToPi;#Lambda(1405) #it{p}_{T} (GeV/c); Bach #pi #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenPtVsBachPtSigmaPlusToPr", "h2GenPtVsBachPtSigmaPlusToPr;#Lambda(1405) #it{p}_{T} (GeV/c); Bach #pi #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenPtVsSigmaMinusPt", "h2GenPtVsBachPtSigmaMinus;#Lambda(1405) #it{p}_{T} (GeV/c); #Sigma #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenPtVsSigmaPlusToPiPt", "h2GenPtVsBachPtSigmaPlusToPi;#Lambda(1405) #it{p}_{T} (GeV/c); #Sigma #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenPtVsSigmaPlusToPrPt", "h2GenPtVsBachPtSigmaPlusToPr;#Lambda(1405) #it{p}_{T} (GeV/c); #Sigma #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenSigmaPtVsKinkPtSigmaMinus", "h2GenSigmaPtVsKinkPtSigmaMinus;#Sigma #it{p}_{T} (GeV/c); Kink #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenSigmaPtVsPiKinkPt", "h2GenSigmaPtVsPiKinkPt;#Sigma #it{p}_{T} (GeV/c); Kink #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2GenSigmaPtVsPrKinkPt", "h2GenSigmaPtVsPrKinkPt;#Sigma #it{p}_{T} (GeV/c); Kink #it{p}_{T} (GeV/c)", {HistType::kTH2F, {ptAxis, ptAxis}});
      rLambda1405.add("h2MassResolutionFromSigmaMinus", "h2MassResolutionFromSigmaMinus", {HistType::kTH2F, {lambda1405MassAxis, massResolutionAxis}});
      rLambda1405.add("h2PtResolutionFromSigmaMinus", "h2PtResolutionFromSigmaMinus", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      rLambda1405.add("h2MassResolutionFromSigmaPlus", "h2MassResolutionFromSigmaPlus", {HistType::kTH2F, {lambda1405MassAxis, massResolutionAxis}});
      rLambda1405.add("h2PtResolutionFromSigmaPlus", "h2PtResolutionFromSigmaPlus", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      rLambda1405.add("h2PtMassMC", "h2PtMassMC", {HistType::kTH2F, {ptAxis, lambda1405MassAxis}});
    }

    // Functional selections
    funcMinQtAlphaAP = TF1("funcMinQtAlphaAP", Form("%s", cutSigmaQtAPMin.value.data()), -1, 1);
    LOGF(info, "funcMinQtAlphaAP: %s", Form("%s", cutSigmaQtAPMin.value.data()));
    funcMaxQtAlphaAP = TF1("funcMaxQtAlphaAP", Form("%s", cutSigmaQtAPMax.value.data()), -1, 1);
    LOGF(info, "funcMaxQtAlphaAP: %s", Form("%s", cutSigmaQtAPMax.value.data()));
    funcMinBachPiPtVsL1405Pt = TF1("funcMinBachPiPtVsL1405Pt", Form("%s", cutMinBachPiPtVsL1405Pt.value.data()), 0., 100);
    LOGF(info, "funcMinBachPiPtVsL1405Pt: %s", Form("%s", cutMinBachPiPtVsL1405Pt.value.data()));
    funcMaxBachPiPtVsL1405Pt = TF1("funcMaxBachPiPtVsL1405Pt", Form("%s", cutMaxBachPiPtVsL1405Pt.value.data()), 0., 100);
    LOGF(info, "funcMaxBachPiPtVsL1405Pt: %s", Form("%s", cutMaxBachPiPtVsL1405Pt.value.data()));
    funcMinSigmaPtVsL1405Pt = TF1("funcMinSigmaPtVsL1405Pt", Form("%s", cutMinSigmaPtVsL1405Pt.value.data()), 0., 100);
    LOGF(info, "funcMinSigmaPtVsL1405Pt: %s", Form("%s", cutMinSigmaPtVsL1405Pt.value.data()));
    funcMaxSigmaPtVsL1405Pt = TF1("funcMaxSigmaPtVsL1405Pt", Form("%s", cutMaxSigmaPtVsL1405Pt.value.data()), 0., 100);
    LOGF(info, "funcMaxSigmaPtVsL1405Pt: %s", Form("%s", cutMaxSigmaPtVsL1405Pt.value.data()));

    rSelections.print();
    rSigmaMinus.print();
    rSigmaPlus.print();
    rLambda1405.print();
  }

  float alphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  float qtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  template <typename TTrack>
  bool selectPiBach(const TTrack& candidate)
  {
    if (std::abs(candidate.tpcNSigmaPi()) > cutNSigTpc) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsBachPi"), 2); // Nsigma Tpc

    if (candidate.tpcNClsFound() < cutNTpcClusPi) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsBachPi"), 3); // Tpc clusters

    if (std::abs(candidate.eta()) > cutEtaDaughter) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsBachPi"), 4); // Eta selection

    return true;
  }

  template <typename TTrack>
  bool selectPiKink(const TTrack& candidate)
  {
    rSelections.fill(HIST("hSelectionsKinkPi"), 0); // All pion kink candidates

    if (std::abs(candidate.tpcNSigmaPi()) > cutNSigTpc) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 1); // Nsigma Tpc

    if (candidate.tpcNClsFound() < cutNTpcClusPi) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 2); // Tpc clusters

    if (std::abs(candidate.eta()) > cutEtaDaughter) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 3); // Eta selection

    if (candidate.itsNCls() < cutNITSClusKink) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 4); // ITS clusters

    if (useTof && !candidate.hasTOF()) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 5); // has Tof

    if (useTof && std::abs(candidate.tofNSigmaPi()) > cutNSigTof) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPi"), 6); // Nsigma Tof

    return true; // Track is selected
  }

  template <typename TTrack>
  bool selectPrKink(const TTrack& candidate)
  {
    rSelections.fill(HIST("hSelectionsKinkPr"), 0); // All proton kink candidates

    if (std::abs(candidate.tpcNSigmaPr()) > cutNSigTpc) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 1); // Nsigma Tpc

    if (candidate.tpcNClsFound() < cutNTpcClusPi) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 2); // Tpc clusters

    if (std::abs(candidate.eta()) > cutEtaDaughter) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 3); // eta selection

    if (candidate.itsNCls() < cutNITSClusKink) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 4); // ITS clusters

    if (useTof && !candidate.hasTOF()) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 5); // has Tof

    if (useTof && std::abs(candidate.tofNSigmaPr()) > cutNSigTof) {
      return false;
    }
    rSelections.fill(HIST("hSelectionsKinkPr"), 6); // Nsigma Tof

    return true; // Track is selected
  }

  template <typename TCand, typename TTrack>
  void fillHistosSigma(const lambda1405candidate& lambda1405Cand, const TCand& sigmaCand, const TTrack& kinkDauTrack)
  {

    if (lambda1405Cand.isSigmaPlus) {
      rSigmaPlus.fill(HIST("h2dPtMassSigmaPlus"), sigmaCand.ptMoth(), sigmaCand.mSigmaPlus());
      rSigmaPlus.fill(HIST("h2dPvContribMassSigmaPlus"), lambda1405Cand.pvContrib, sigmaCand.mSigmaPlus());
      rSigmaPlus.fill(HIST("h2dPvContribPtSigmaPlus"), lambda1405Cand.pvContrib, sigmaCand.ptMoth());
      rSigmaPlus.fill(HIST("hMassXiMinusSigmaPlus"), sigmaCand.mXiMinus(), sigmaCand.mSigmaPlus());
      rSigmaPlus.fill(HIST("hSigmaPlusArmPod"), lambda1405Cand.sigmaAlphaAP, lambda1405Cand.sigmaQtAP);
      rSigmaPlus.fill(HIST("hSigmaPlusRadius"), lambda1405Cand.sigmaRadius, sigmaCand.ptMoth());
      rSigmaPlus.fill(HIST("hSigmaPlusDcaToPv"), sigmaCand.dcaMothPv(), sigmaCand.ptMoth());
      rSigmaPlus.fill(HIST("hSigmaPlusDcaKinkDauToPv"), sigmaCand.dcaDaugPv(), sigmaCand.ptMoth());
      rSigmaPlus.fill(HIST("h2KinkPrPtNSigTofSigmaPlus"), lambda1405Cand.kinkPt, lambda1405Cand.kinkPrNSigTof);
      // Fill QA histos for kink daughter
      rSigmaPlus.fill(HIST("hSigmaPlusKinkTpcNSigPr"), kinkDauTrack.tpcNSigmaPr(), kinkDauTrack.pt());
      rSigmaPlus.fill(HIST("hSigmaPlusKinkTofNSigPr"), kinkDauTrack.tofNSigmaPr(), kinkDauTrack.pt());
    } else {
      rSigmaMinus.fill(HIST("h2dPtMassSigmaMinus"), sigmaCand.ptMoth(), sigmaCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2dPvContribMassSigmaMinus"), lambda1405Cand.pvContrib, sigmaCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2dPvContribPtSigmaMinus"), lambda1405Cand.pvContrib, sigmaCand.ptMoth());
      rSigmaMinus.fill(HIST("hMassXiMinusSigmaMinus"), sigmaCand.mXiMinus(), sigmaCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("hSigmaMinusArmPod"), lambda1405Cand.sigmaAlphaAP, lambda1405Cand.sigmaQtAP);
      rSigmaMinus.fill(HIST("hSigmaMinusRadius"), lambda1405Cand.sigmaRadius, sigmaCand.ptMoth());
      rSigmaMinus.fill(HIST("hSigmaMinusDcaToPv"), sigmaCand.dcaMothPv(), sigmaCand.ptMoth());
      rSigmaMinus.fill(HIST("hSigmaMinusDcaKinkDauToPv"), sigmaCand.dcaDaugPv(), sigmaCand.ptMoth());
      rSigmaMinus.fill(HIST("h2KinkPiPtNSigTofSigmaMinus"), lambda1405Cand.kinkPt, lambda1405Cand.kinkPiNSigTof);
      // Fill QA histos for kink daughter
      rSigmaMinus.fill(HIST("hSigmaMinusKinkTpcNSigPi"), kinkDauTrack.tpcNSigmaPi(), kinkDauTrack.pt());
      rSigmaMinus.fill(HIST("hSigmaMinusKinkTofNSigPi"), kinkDauTrack.tofNSigmaPi(), kinkDauTrack.pt());
    }
  }

  template <bool IsMC>
  void fillHistosLambda1405(const lambda1405candidate& cand)
  {

    // Fill QA histos for Lambda(1405) candidate
    rLambda1405.fill(HIST("hMassL1405"), cand.massL1405, cand.pt());
    rLambda1405.fill(HIST("hMassXi1530"), cand.massXi1530, cand.pt());
    rLambda1405.fill(HIST("hPx"), cand.px);
    rLambda1405.fill(HIST("hPy"), cand.py);
    rLambda1405.fill(HIST("hPz"), cand.pz);
    rLambda1405.fill(HIST("hPt"), cand.pt());
    rLambda1405.fill(HIST("hPhi"), cand.phi);

    // Bachelor Pi
    rLambda1405.fill(HIST("h2BachPiPtNSigTpc"), cand.bachPiPt, cand.bachPiNSigTpc);
    rLambda1405.fill(HIST("h2BachPiPtNSigTof"), cand.bachPiPt, cand.bachPiNSigTof);

    // Fill sparse histo
    std::vector<double> sparseEntry = {cand.massL1405, cand.pt(), cand.sigmaMinusMass, cand.sigmaAlphaAP, cand.sigmaQtAP, cand.dcaSigmaToPv, cand.kinkDcaDauToPv};
    if constexpr (IsMC) {
      sparseEntry.push_back(cand.pvContrib);
    }
    if (doprocessDataWCentQVecs) {
      sparseEntry.push_back(cand.cent);
      sparseEntry.push_back(cand.scalarProd);
    }
    auto hSparse = rLambda1405.get<THnSparse>(HIST("hSparseL1405"));
    hSparse->Fill(sparseEntry.data());
  }

  template <typename TColl>
  void constructCollCandidates(const TColl& collision, aod::KinkCands::iterator const& sigmaCand, TracksFull const& tracks, std::vector<lambda1405candidate>& selectedCandidates)
  {

    lambda1405candidate lambda1405Cand{};

    rSelections.fill(HIST("hSelectionsL1405"), 0);      // All candidates
    rSelections.fill(HIST("hSelectionsSigmaMinus"), 0); // All Sigma- candidates
    rSelections.fill(HIST("hSelectionsSigmaPlus"), 0);  // All Sigma+ candidates

    auto kinkDauTrack = sigmaCand.template trackDaug_as<TracksFull>();
    bool isPiKink = selectPiKink(kinkDauTrack);
    bool isPrKink = selectPrKink(kinkDauTrack);
    if (!isPiKink && !isPrKink) {
      return;
    }

    if (isPiKink) { // Dominated by Sigma-, Sigma+ treated as contamination
      lambda1405Cand.isSigmaMinus = true;
      lambda1405Cand.isSigmaPlus = false;
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 1); // Passed kink sel
    } else {                                              // Only Sigma+ can have a proton as kink daughter
      lambda1405Cand.isSigmaMinus = false;
      lambda1405Cand.isSigmaPlus = true;
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 1); // Passed kink sel
    }

    // Sigma- or AntiSigma+ candidates
    if (lambda1405Cand.isSigmaMinus) {
      rSigmaMinus.fill(HIST("h2PtMassSigmaMinusBeforeCuts"), sigmaCand.ptMoth(), sigmaCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2PtPiKinkNSigBeforeCutsSigmaMinus"), kinkDauTrack.pt(), kinkDauTrack.tpcNSigmaPi());
    }
    if (lambda1405Cand.isSigmaPlus) {
      rSigmaPlus.fill(HIST("h2PtMassSigmaPlusBeforeCuts"), sigmaCand.ptMoth(), sigmaCand.mSigmaPlus());
      rSigmaPlus.fill(HIST("h2PtPrKinkNSigBeforeCuts"), kinkDauTrack.pt(), kinkDauTrack.tpcNSigmaPr());
    }

    if (lambda1405Cand.isSigmaMinus && (sigmaCand.mSigmaMinus() < o2::constants::physics::MassSigmaMinus - cutSigmaMass ||
                                        sigmaCand.mSigmaMinus() > o2::constants::physics::MassSigmaMinus + cutSigmaMass)) {
      return;
    }
    if (lambda1405Cand.isSigmaPlus && (sigmaCand.mSigmaPlus() < o2::constants::physics::MassSigmaPlus - cutSigmaMass ||
                                       sigmaCand.mSigmaPlus() > o2::constants::physics::MassSigmaPlus + cutSigmaMass)) {
      return;
    }
    if (lambda1405Cand.isSigmaMinus) {
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 2); // Passed mass sel
    } else {
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 2); // Passed mass sel
    }

    if (std::abs(sigmaCand.dcaMothPv()) > cutDCAtoPvSigma) {
      return;
    }
    if (lambda1405Cand.isSigmaMinus) {
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 3); // Passed cutDCAtoPvSigma
    } else {
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 3); // Passed cutDCAtoPvSigma
    }

    if (std::abs(sigmaCand.dcaDaugPv()) < cutDCAtoPvPiFromSigma) {
      return;
    }
    if (lambda1405Cand.isSigmaMinus) {
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 4); // cutDCAtoPvPiFromSigma
    } else {
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 4); // cutDCAtoPvPiFromSigma
    }

    float sigmaRad = std::hypot(sigmaCand.xDecVtx(), sigmaCand.yDecVtx());
    if (sigmaRad < cutSigmaRadius) {
      return;
    }
    if (lambda1405Cand.isSigmaMinus) {
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 5); // Passed radius sel
    } else {
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 5); // Passed radius sel
    }

    auto kinkDauMom = std::array{sigmaCand.pxDaug(), sigmaCand.pyDaug(), sigmaCand.pzDaug()};
    auto sigmaMom = std::array{sigmaCand.pxMoth(), sigmaCand.pyMoth(), sigmaCand.pzMoth()};
    // Sigma properties
    lambda1405Cand.sigmaId = sigmaCand.globalIndex();
    lambda1405Cand.sigmaMinusMass = sigmaCand.mSigmaMinus();
    lambda1405Cand.sigmaPlusMass = sigmaCand.mSigmaPlus();
    lambda1405Cand.xiMinusMass = sigmaCand.mXiMinus();
    lambda1405Cand.sigmaSign = sigmaCand.mothSign();
    lambda1405Cand.sigmaAlphaAP = alphaAP(sigmaMom, kinkDauMom);
    lambda1405Cand.sigmaQtAP = qtAP(sigmaMom, kinkDauMom);
    lambda1405Cand.sigmaPt = sigmaCand.ptMoth();
    lambda1405Cand.sigmaRadius = sigmaRad;
    lambda1405Cand.dcaSigmaToPv = sigmaCand.dcaMothPv();

    if (lambda1405Cand.sigmaQtAP < funcMinQtAlphaAP.Eval(lambda1405Cand.sigmaAlphaAP) ||
        lambda1405Cand.sigmaQtAP > funcMaxQtAlphaAP.Eval(lambda1405Cand.sigmaAlphaAP)) {
      return;
    }
    if (lambda1405Cand.isSigmaMinus) {
      rSelections.fill(HIST("hSelectionsSigmaMinus"), 6); // Passed AP sel
    } else {
      rSelections.fill(HIST("hSelectionsSigmaPlus"), 6); // Passed AP sel
    }

    // Kink daughter properties
    lambda1405Cand.kinkDauId = kinkDauTrack.globalIndex();
    lambda1405Cand.kinkPt = kinkDauTrack.pt();
    lambda1405Cand.kinkPiNSigTpc = kinkDauTrack.tpcNSigmaPi();
    lambda1405Cand.kinkPiNSigTof = kinkDauTrack.tofNSigmaPi();
    lambda1405Cand.kinkPrNSigTpc = kinkDauTrack.tpcNSigmaPr();
    lambda1405Cand.kinkPrNSigTof = kinkDauTrack.tofNSigmaPr();
    lambda1405Cand.kinkDcaDauToPv = sigmaCand.dcaDaugPv();

    rSelections.fill(HIST("hSelectionsL1405"), 1); // Passed Sigma sel

    // Collision properties
    lambda1405Cand.pvContrib = collision.numContrib();
    if constexpr (requires { collision.centFT0C(); }) {
      lambda1405Cand.cent = collision.centFT0C();
    } else {
      lambda1405Cand.cent = -1; // Not available
    }

    fillHistosSigma(lambda1405Cand, sigmaCand, kinkDauTrack);

    for (const auto& piTrack : tracks) {
      rSelections.fill(HIST("hSelectionsBachPi"), 0); // All bachelors

      // Needed to avoid spurious correlations in the Like-Sign case
      if (piTrack.globalIndex() == kinkDauTrack.globalIndex()) {
        continue; // Skip the kink daughter track
      }

      bool acceptPair{false};
      if (doLikeSignBkg) {
        acceptPair = (piTrack.sign() == sigmaCand.mothSign());
      } else {
        acceptPair = (piTrack.sign() != sigmaCand.mothSign());
      }
      if (!acceptPair) {
        continue;
      }
      rSelections.fill(HIST("hSelectionsBachPi"), 1);

      if (!selectPiBach(piTrack)) {
        continue;
      }
      rSelections.fill(HIST("hSelectionsBachPi"), 5); // PID sel
      rSelections.fill(HIST("hSelectionsL1405"), 2);  // Bach Pi selection

      auto piMom = std::array{piTrack.px(), piTrack.py(), piTrack.pz()};
      float invMass{-1.f};
      if (lambda1405Cand.isSigmaMinus) {
        invMass = RecoDecay::m(std::array{sigmaMom, piMom}, std::array{o2::constants::physics::MassSigmaMinus, o2::constants::physics::MassPiPlus});
      } else if (lambda1405Cand.isSigmaPlus) {
        invMass = RecoDecay::m(std::array{sigmaMom, piMom}, std::array{o2::constants::physics::MassSigmaPlus, o2::constants::physics::MassPiMinus});
      }
      if (invMass > cutUpperMass) {
        continue;
      }
      rSelections.fill(HIST("hSelectionsL1405"), 3); // Upper mass selection

      // Daughter Pi properties
      lambda1405Cand.bachPiId = piTrack.globalIndex();
      lambda1405Cand.bachPiPt = piTrack.pt();
      lambda1405Cand.bachPiNSigTpc = piTrack.tpcNSigmaPi();
      if (useTof) {
        lambda1405Cand.bachPiNSigTof = piTrack.tofNSigmaPi();
      } else {
        lambda1405Cand.bachPiNSigTof = -999; // Not used if Tof is not enabled
      }

      // Lambda(1405) candidate properties
      lambda1405Cand.massL1405 = invMass;
      lambda1405Cand.massXi1530 = RecoDecay::m(std::array{sigmaMom, piMom}, std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
      lambda1405Cand.px = sigmaMom[0] + piMom[0];
      lambda1405Cand.py = sigmaMom[1] + piMom[1];
      lambda1405Cand.pz = sigmaMom[2] + piMom[2];

      lambda1405Cand.phi = std::atan2(lambda1405Cand.py, lambda1405Cand.px);
      lambda1405Cand.scalarProd = -1;
      if constexpr (requires { collision.qvecFT0CRe(); }) {
        float const xQVec = collision.qvecFT0CRe();
        float const yQVec = collision.qvecFT0CIm();
        float const cos2Phi = std::cos(2 * lambda1405Cand.phi);
        float const sin2Phi = std::sin(2 * lambda1405Cand.phi);
        lambda1405Cand.scalarProd = cos2Phi * xQVec + sin2Phi * yQVec;
      }

      // Check correlations between transverse momenta of L1405, sigma and bachelor pi
      if (std::hypot(sigmaCand.pxMoth(), sigmaCand.pyMoth()) < funcMinSigmaPtVsL1405Pt.Eval(lambda1405Cand.pt()) ||
          std::hypot(sigmaCand.pxMoth(), sigmaCand.pyMoth()) > funcMaxSigmaPtVsL1405Pt.Eval(lambda1405Cand.pt())) {
        continue;
      }

      if (piTrack.pt() < funcMinBachPiPtVsL1405Pt.Eval(lambda1405Cand.pt()) ||
          piTrack.pt() > funcMaxBachPiPtVsL1405Pt.Eval(lambda1405Cand.pt())) {
        continue;
      }
      rSelections.fill(HIST("hSelectionsL1405"), 4); // Pt correlations

      if (lambda1405Cand.pt() < cutMinPtL1405) {
        continue;
      }
      rSelections.fill(HIST("hSelectionsL1405"), 5); // Accepted

      selectedCandidates.push_back(lambda1405Cand);
    }
  }

  template <typename mcTrack>
  bool checkSigmaKinkMC(const mcTrack& mcTrackSigma, const mcTrack& mcTrackKinkDau, float sigmaAbsPDG, float kinkAbsPDG, aod::McParticles const&)
  {
    if (std::abs(mcTrackSigma.pdgCode()) != sigmaAbsPDG || std::abs(mcTrackKinkDau.pdgCode()) != kinkAbsPDG) {
      return false; // Not a valid Sigma kink decay
    }
    if (!mcTrackKinkDau.has_mothers()) {
      return false; // No mothers found
    }
    // Check if the kink comes from the Sigma
    bool isKinkFromSigma = false;
    for (const auto& mcMother : mcTrackKinkDau.template mothers_as<aod::McParticles>()) {
      if (mcMother.globalIndex() == mcTrackSigma.globalIndex()) {
        isKinkFromSigma = true;
        break;
      }
    }
    return isKinkFromSigma; // Return true if the kink comes from the Sigma
  }

  template <typename TCollision, typename TCand, typename TTrack>
  void fillOutputData(const TCollision& collision, const TCand& sigmaCands, const TTrack& tracks)
  {
    if constexpr (requires { collision.centFT0C(); }) {
      if (collision.centFT0C() < centralityMin || collision.centFT0C() > centralityMax) {
        return;
      }
    }
    if (std::abs(collision.posZ()) > cutZVertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& sigmaCand : sigmaCands) {
      std::vector<lambda1405candidate> selectedCandidates;
      constructCollCandidates(collision, sigmaCand, tracks, selectedCandidates);
      for (auto& lambda1405Cand : selectedCandidates) {
        if (lambda1405Cand.isSigmaMinus) {
          rLambda1405.fill(HIST("h2SigmaMinusMassVsLambdaMass"), lambda1405Cand.massL1405, lambda1405Cand.sigmaMinusMass);
        } else {
          rLambda1405.fill(HIST("h2SigmaPlusMassVsLambdaMass"), lambda1405Cand.massL1405, lambda1405Cand.sigmaPlusMass);
        }
        if (fillOutputTree) {
          float const ptCand = lambda1405Cand.pt();
          if (downSampleFactor < 1.) {
            float const pseudoRndm = ptCand * 1000. - static_cast<int64_t>(ptCand * 1000);
            if (ptCand < ptDownSampleMax && pseudoRndm >= downSampleFactor) {
              continue;
            }
          }
          outputDataTable(lambda1405Cand.px, lambda1405Cand.py, lambda1405Cand.pz,
                          lambda1405Cand.massL1405, lambda1405Cand.massXi1530,
                          lambda1405Cand.sigmaMinusMass, lambda1405Cand.sigmaPlusMass, lambda1405Cand.xiMinusMass,
                          lambda1405Cand.sigmaPt, lambda1405Cand.sigmaAlphaAP, lambda1405Cand.sigmaQtAP, lambda1405Cand.sigmaRadius,
                          lambda1405Cand.kinkPt,
                          lambda1405Cand.kinkPiNSigTpc, lambda1405Cand.kinkPiNSigTof,
                          lambda1405Cand.kinkPrNSigTpc, lambda1405Cand.kinkPrNSigTof,
                          lambda1405Cand.kinkDcaDauToPv,
                          lambda1405Cand.bachPiNSigTpc, lambda1405Cand.bachPiNSigTof);
        }
        if constexpr (requires { collision.qvecFT0CRe(); }) {
          float const xQVec = collision.qvecFT0CRe();
          float const yQVec = collision.qvecFT0CIm();
          float const cos2Phi = std::cos(2 * lambda1405Cand.phi);
          float const sin2Phi = std::sin(2 * lambda1405Cand.phi);
          lambda1405Cand.scalarProd = cos2Phi * xQVec + sin2Phi * yQVec;
          float const ptCand = lambda1405Cand.pt();
          rLambda1405.fill(HIST("hScalarProd"), lambda1405Cand.scalarProd, ptCand);
          if (fillFlowTree) {
            if (downSampleFactor < 1.) {
              float const pseudoRndm = ptCand * 1000. - static_cast<int64_t>(ptCand * 1000);
              if (ptCand < ptDownSampleMax && pseudoRndm >= downSampleFactor) {
                continue;
              }
            }
            outputDataFlowTable(ptCand, lambda1405Cand.massL1405,
                                lambda1405Cand.sigmaPt,
                                lambda1405Cand.sigmaMinusMass, lambda1405Cand.sigmaPlusMass,
                                lambda1405Cand.sigmaAlphaAP, lambda1405Cand.sigmaQtAP,
                                lambda1405Cand.kinkPiNSigTpc, lambda1405Cand.kinkPiNSigTof,
                                lambda1405Cand.kinkPrNSigTpc, lambda1405Cand.kinkPrNSigTof,
                                lambda1405Cand.kinkDcaDauToPv,
                                lambda1405Cand.bachPiNSigTpc, lambda1405Cand.bachPiNSigTof,
                                lambda1405Cand.scalarProd, collision.centFT0C());
          }
        }
        fillHistosLambda1405<false>(lambda1405Cand);
      }
    }
  }

  void processData(CollisionsFull::iterator const& collision,
                   aod::KinkCands const& kinkCands,
                   TracksFull const& tracks)
  {
    fillOutputData(collision, kinkCands, tracks);
  }
  PROCESS_SWITCH(lambda1405analysis, processData, "Data processing", true);

  void processDataWCentQVecs(CollisionsCentSel::iterator const& collision,
                             aod::KinkCands const& kinkCands,
                             TracksFull const& tracks)
  {
    fillOutputData(collision, kinkCands, tracks);
  }
  PROCESS_SWITCH(lambda1405analysis, processDataWCentQVecs, "Data processing with centrality and Q vectors info", false);

  template <typename TMother>
  int matchGenDecay(const TMother& motherPart, const aod::McParticles& mcParticles, std::array<int, 3>& daugsIdxs)
  {
    int pdgMother = motherPart.pdgCode();
    int8_t sign = 0;

    int MaxDepthFinState = 2; // Maximum depth to look for the decay chain

    std::vector<int> arrL1405Daugs = {};
    std::array<int, 3> finalState;
    int decayChannel = -1;

    // Match L(1405) --> n pi- pi+ final state
    if (pdgMother > 0) { // Change sign of neutral decay products
      finalState = {PDG_t::kNeutron, PDG_t::kPiMinus, PDG_t::kPiPlus};
    } else {
      finalState = {-PDG_t::kNeutron, PDG_t::kPiMinus, PDG_t::kPiPlus};
    }
    if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, motherPart, lambda1405PdgCode, finalState, true, &sign, MaxDepthFinState)) {
      // Match the intermediate Sigma
      RecoDecay::getDaughters<true>(motherPart, &arrL1405Daugs, std::array{0}, 1);
      for (auto iProng = 0u; iProng < arrL1405Daugs.size(); ++iProng) {
        auto daughI = mcParticles.rawIteratorAt(arrL1405Daugs[iProng]);
        if (std::abs(daughI.pdgCode()) == PDG_t::kSigmaPlus) {
          decayChannel = kSigmaPlusPiToPiPiNeutron;
          daugsIdxs[1] = arrL1405Daugs[iProng];
        }
        if (std::abs(daughI.pdgCode()) == PDG_t::kSigmaMinus) {
          decayChannel = kSigmaMinusPiToPiPiNeutron;
          daugsIdxs[1] = arrL1405Daugs[iProng];
        }
        if (std::abs(daughI.pdgCode()) == PDG_t::kPiPlus) {
          daugsIdxs[0] = arrL1405Daugs[iProng];
        }
      }
    }

    // Match L(1405) --> p pi0 pi+ final state, only possible for Sigma+
    if (pdgMother > 0) { // Change sign of neutral decay products
      finalState = {PDG_t::kProton, PDG_t::kPi0, PDG_t::kPiMinus};
    } else {
      finalState = {PDG_t::kProton, -PDG_t::kPi0, PDG_t::kPiMinus};
    }

    if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, motherPart, lambda1405PdgCode, finalState, true, &sign, MaxDepthFinState)) {
      // Match the intermediate Sigma
      RecoDecay::getDaughters<true>(motherPart, &arrL1405Daugs, std::array{0}, 1);
      for (auto iProng = 0u; iProng < arrL1405Daugs.size(); ++iProng) {
        auto daughI = mcParticles.rawIteratorAt(arrL1405Daugs[iProng]);
        if (std::abs(daughI.pdgCode()) == PDG_t::kSigmaPlus) {
          decayChannel = kSigmaPlusPiToPiPiProton;
          daugsIdxs[1] = arrL1405Daugs[iProng];
        }
        if (std::abs(daughI.pdgCode()) == PDG_t::kPiPlus) {
          daugsIdxs[0] = arrL1405Daugs[iProng];
        }
      }
    }

    std::vector<int> arrSigmaDaugs = {};
    auto sigmaDaug = mcParticles.rawIteratorAt(daugsIdxs[1]);
    RecoDecay::getDaughters<true>(sigmaDaug, &arrSigmaDaugs, std::array{0}, 1);
    for (auto iProng = 0u; iProng < arrSigmaDaugs.size(); ++iProng) {
      auto daughSigma = mcParticles.rawIteratorAt(arrSigmaDaugs[iProng]);
      if (std::abs(daughSigma.pdgCode()) == PDG_t::kPiPlus || std::abs(daughSigma.pdgCode()) == PDG_t::kProton) {
        daugsIdxs[2] = arrSigmaDaugs[iProng];
      }
    }

    return decayChannel;
  }

  template <typename TCollision, typename TTrack>
  void fillOutputMc(const TCollision& recoCollisions,
                    const aod::KinkCands& sigmaCands,
                    const aod::McTrackLabels& trackLabelsMC,
                    const TTrack& tracks,
                    const aod::McParticles& particlesMC)
  {
    for (const auto& collision : recoCollisions) {
      if constexpr (requires { collision.centFT0C(); }) {
        auto cent = collision.centFT0C();
        if (cent < centralityMin || cent > centralityMax) {
          return;
        }
      }
      if (std::abs(collision.posZ()) > cutZVertex) { // || !collision.sel8()) {
        continue;
      }
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto sigmaCandsPerCol = sigmaCands.sliceBy(mKinkPerCol, collision.globalIndex());
      auto tracksPerCol = tracks.sliceBy(mPerColTracks, collision.globalIndex());
      for (const auto& sigmaCand : sigmaCandsPerCol) {

        // Perform the sigma matching here, so sigma histograms
        // can be used for efficiency studies on the sigma
        auto labelSigma = trackLabelsMC.rawIteratorAt(sigmaCand.trackMothId());
        auto labelKinkDaug = trackLabelsMC.rawIteratorAt(sigmaCand.trackDaugId());
        if (!labelSigma.has_mcParticle() || !labelKinkDaug.has_mcParticle()) {
          continue; // Skip if no valid MC association
        }
        auto genSigma = labelSigma.template mcParticle_as<aod::McParticles>();
        auto genKinkDaug = labelKinkDaug.template mcParticle_as<aod::McParticles>();

        bool isSigmaMinusKink = checkSigmaKinkMC(genSigma, genKinkDaug, sigmaMinusPdgCode, 211, particlesMC);
        bool isSigmaPlusToPiKink = checkSigmaKinkMC(genSigma, genKinkDaug, sigmaPlusPdgCode, 211, particlesMC);
        bool isSigmaPlusToPrKink = checkSigmaKinkMC(genSigma, genKinkDaug, sigmaPlusPdgCode, 2212, particlesMC);

        if (!isSigmaMinusKink && !isSigmaPlusToPiKink && !isSigmaPlusToPrKink) {
          continue; // Skip if not a valid Sigma kink decay
        }

        std::vector<lambda1405candidate> selectedCandidates;
        constructCollCandidates(collision, sigmaCand, tracksPerCol, selectedCandidates);
        for (auto& lambda1405Cand : selectedCandidates) {
          rLambda1405.fill(HIST("hRecoL1405"), 0., lambda1405Cand.pt()); // All reconstructed

          // Do MC association
          auto labelBachPi = trackLabelsMC.rawIteratorAt(lambda1405Cand.bachPiId);
          if (!labelBachPi.has_mcParticle()) {
            continue; // Skip if no valid MC association
          }

          auto genBachPi = labelBachPi.template mcParticle_as<aod::McParticles>();
          if (std::abs(genBachPi.pdgCode()) != PDG_t::kPiPlus) {
            continue; // Skip if not a valid pion candidate
          }
          rLambda1405.fill(HIST("hRecoL1405"), 1., lambda1405Cand.pt()); // Has bach pi

          if (!genSigma.has_mothers() || !genBachPi.has_mothers()) {
            continue; // Skip if no mothers found
          }
          rLambda1405.fill(HIST("hRecoL1405"), 2., lambda1405Cand.pt()); // Has mothers for Sigma and Pi

          // check that labpi and labsigma have the same mother (a lambda1405 candidate)
          int lambda1405Id = -1;
          for (const auto& piMother : genBachPi.template mothers_as<aod::McParticles>()) {
            for (const auto& sigmaMother : genSigma.template mothers_as<aod::McParticles>()) {
              if (piMother.globalIndex() == sigmaMother.globalIndex() && std::abs(piMother.pdgCode()) == lambda1405PdgCode) {
                lambda1405Id = piMother.globalIndex();
                break; // Found the mother, exit loop
              }
            }
          }
          if (lambda1405Id == -1) {
            continue; // Skip if the Sigma and pion do not share the same lambda1405 candidate
          }
          rLambda1405.fill(HIST("hRecoL1405"), 3., lambda1405Cand.pt()); // Has same mother

          auto lambda1405Mother = particlesMC.rawIteratorAt(lambda1405Id);
          float lambda1405Mass = std::sqrt(lambda1405Mother.e() * lambda1405Mother.e() - lambda1405Mother.p() * lambda1405Mother.p());
          if (lambda1405Cand.isSigmaMinus) {
            rLambda1405.fill(HIST("h2SigmaMinusMassVsLambdaMass"), lambda1405Cand.massL1405, lambda1405Cand.sigmaMinusMass);
            rLambda1405.fill(HIST("h2MassResolutionFromSigmaMinus"), lambda1405Mass, lambda1405Mass - lambda1405Cand.massL1405);
            rLambda1405.fill(HIST("h2PtResolutionFromSigmaMinus"), lambda1405Cand.pt(), lambda1405Cand.pt() - lambda1405Mother.pt());
          }
          if (lambda1405Cand.isSigmaPlus) {
            rLambda1405.fill(HIST("h2SigmaPlusMassVsLambdaMass"), lambda1405Cand.massL1405, lambda1405Cand.sigmaPlusMass);
            rLambda1405.fill(HIST("h2MassResolutionFromSigmaPlus"), lambda1405Mass, lambda1405Mass - lambda1405Cand.massL1405);
            rLambda1405.fill(HIST("h2PtResolutionFromSigmaPlus"), lambda1405Cand.pt(), lambda1405Cand.pt() - lambda1405Mother.pt());
          }

          if (fillOutputTree) {
            outputDataTableMC(lambda1405Cand.px, lambda1405Cand.py, lambda1405Cand.pz,
                              lambda1405Cand.massL1405, lambda1405Cand.massXi1530,
                              lambda1405Cand.sigmaMinusMass, lambda1405Cand.sigmaPlusMass, lambda1405Cand.xiMinusMass,
                              lambda1405Cand.sigmaPt, lambda1405Cand.sigmaAlphaAP, lambda1405Cand.sigmaQtAP, lambda1405Cand.sigmaRadius,
                              lambda1405Cand.kinkPt,
                              lambda1405Cand.kinkPiNSigTpc, lambda1405Cand.kinkPiNSigTof,
                              lambda1405Cand.kinkPrNSigTpc, lambda1405Cand.kinkPrNSigTof,
                              lambda1405Cand.kinkDcaDauToPv,
                              lambda1405Cand.bachPiNSigTpc, lambda1405Cand.bachPiNSigTof,
                              lambda1405Mother.pt(), lambda1405Mass, genSigma.pdgCode(), genBachPi.pdgCode());
          }
          fillHistosLambda1405<true>(lambda1405Cand);
        }
      }
    }

    // Loop over generated particles to fill MC histograms
    for (const auto& mcPart : particlesMC) {
      if (std::abs(mcPart.pdgCode()) != lambda1405PdgCode &&
          std::abs(mcPart.pdgCode()) != sigmaMinusPdgCode &&
          std::abs(mcPart.pdgCode()) != sigmaPlusPdgCode) {
        continue; // Only consider Lambda(1405) and Sigma candidates
      }

      // Compute generated PV contributors
      const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, mcPart.mcCollision().globalIndex());
      if (recoCollsPerMcColl.size() == 0) {
        continue; // Skip if no reconstructed collisions associated with this MC collision
      }
      unsigned maxNumContrib = 0;
      for (const auto& recCol : recoCollsPerMcColl) {
        maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
      }

      // Needed for Sigma efficiency vs PV contributors
      if (std::abs(mcPart.pdgCode()) == sigmaMinusPdgCode) {
        rSigmaMinus.fill(HIST("h2GenSigmaMinusPvContribPt"), maxNumContrib, mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == sigmaPlusPdgCode) {
        rSigmaPlus.fill(HIST("h2GenSigmaPlusPvContribPt"), maxNumContrib, mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) != lambda1405PdgCode) {
        continue; // Only consider Lambda(1405) candidates
      }
      std::array<int, 3> idxsSigmaKinkBachPi{-1, -1, -1};
      int decayChannel = matchGenDecay(mcPart, particlesMC, idxsSigmaKinkBachPi);
      if (decayChannel == -1) {
        continue; // Skip if it doesn't match the decay channels of interest
      }
      rLambda1405.fill(HIST("hGenL1405"), decayChannel, mcPart.pt());

      auto bachPi = particlesMC.rawIteratorAt(idxsSigmaKinkBachPi[0]);
      auto sigmaDaug = particlesMC.rawIteratorAt(idxsSigmaKinkBachPi[1]);
      auto sigmaKinkDaug = particlesMC.rawIteratorAt(idxsSigmaKinkBachPi[2]);
      // Generated Armenteros-Podolanski variables
      float genSigmaAlphaAP = alphaAP({sigmaDaug.px(), sigmaDaug.py(), sigmaDaug.pz()}, {sigmaKinkDaug.px(), sigmaKinkDaug.py(), sigmaKinkDaug.pz()});
      float genSigmaQtAP = qtAP({sigmaDaug.px(), sigmaDaug.py(), sigmaDaug.pz()}, {sigmaKinkDaug.px(), sigmaKinkDaug.py(), sigmaKinkDaug.pz()});
      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      rLambda1405.fill(HIST("h2PtMassMC"), mcPart.pt(), mcMass);
      rLambda1405.fill(HIST("h2GenL1405PvContribPt"), maxNumContrib, mcPart.pt());
      if (decayChannel == kSigmaMinusPiToPiPiNeutron) {
        rLambda1405.fill(HIST("h2GenSigmaMinusArmPod"), genSigmaAlphaAP, genSigmaQtAP);
        rLambda1405.fill(HIST("h2GenPtVsBachPtSigmaMinus"), mcPart.pt(), bachPi.pt());
        rLambda1405.fill(HIST("h2GenPtVsSigmaMinusPt"), mcPart.pt(), sigmaDaug.pt());
        rLambda1405.fill(HIST("h2GenSigmaPtVsKinkPtSigmaMinus"), sigmaDaug.pt(), sigmaKinkDaug.pt());
      }
      if (decayChannel == kSigmaPlusPiToPiPiNeutron) {
        rLambda1405.fill(HIST("h2GenSigmaPlusArmPod"), genSigmaAlphaAP, genSigmaQtAP);
        rLambda1405.fill(HIST("h2GenPtVsBachPtSigmaPlusToPi"), mcPart.pt(), bachPi.pt());
        rLambda1405.fill(HIST("h2GenPtVsSigmaPlusToPiPt"), mcPart.pt(), sigmaDaug.pt());
        rLambda1405.fill(HIST("h2GenSigmaPtVsPiKinkPt"), sigmaDaug.pt(), sigmaKinkDaug.pt());
      }
      if (decayChannel == kSigmaPlusPiToPiPiProton) {
        rLambda1405.fill(HIST("h2GenSigmaPlusArmPod"), genSigmaAlphaAP, genSigmaQtAP);
        rLambda1405.fill(HIST("h2GenPtVsBachPtSigmaPlusToPr"), mcPart.pt(), bachPi.pt());
        rLambda1405.fill(HIST("h2GenPtVsSigmaPlusToPrPt"), mcPart.pt(), sigmaDaug.pt());
        rLambda1405.fill(HIST("h2GenSigmaPtVsPrKinkPt"), sigmaDaug.pt(), sigmaKinkDaug.pt());
      }
    }
  }

  void processMc(CollisionsFullMc const& recoCollisions,
                 aod::KinkCands const& kinkCands,
                 aod::McTrackLabels const& trackLabelsMC,
                 aod::McParticles const& particlesMC,
                 TracksFull const& tracks,
                 aod::McCollisions const&)
  {
    fillOutputMc(recoCollisions, kinkCands, trackLabelsMC, tracks, particlesMC);
  }
  PROCESS_SWITCH(lambda1405analysis, processMc, "MC processing", false);

  void processMcWCentSel(McRecoCollisionsCentSel const& recoCollisions,
                         aod::KinkCands const& kinkCands,
                         aod::McTrackLabels const& trackLabelsMC,
                         aod::McParticles const& particlesMC,
                         TracksFull const& tracks)
  {
    fillOutputMc(recoCollisions, kinkCands, trackLabelsMC, tracks, particlesMC);
  }
  PROCESS_SWITCH(lambda1405analysis, processMcWCentSel, "MC processing with centrality selection", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambda1405analysis>(cfgc)};
}
