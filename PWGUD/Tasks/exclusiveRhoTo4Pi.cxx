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

/// \file exclusiveRhoTo4Pi.cxx
/// \brief Task for analyzing exclusive rho decays to 4 pions
/// \author Anantha Padmanabhan M Nair

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include <TMath.h>
#include <TString.h>

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using PtEtaPhiMVector = ROOT::Math::PtEtaPhiMVector;
using Boost = ROOT::Math::Boost;
using XYZVectorF = ROOT::Math::XYZVectorF;
using PxPyPzEVector = ROOT::Math::PxPyPzEVector;
using PxPyPzMVector = ROOT::Math::PxPyPzMVector;

struct ExclusiveRhoTo4Pi {
  SGSelector sgSelector;
  // Numbers for background estimation
  int zero = 0;
  int one = 1;
  int two = 2;
  int three = 3;
  int four = 4;
  // PDG Codes and rho mass
  double mRho0 = 0.77526; // GeV/c^2
  int rhoPrime = 30113;
  // Pb-Pb at 5.36 TeV
  double halfSqrtSnn = 2680.;
  double massOfLead208 = 193.6823;
  double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);
  // Run Numbers
  static int runNos[113];
  static int numRunNums;
  // Histogram Registry
  HistogramRegistry histosDataCounter{"Counters", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosQA{"QA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosPID{"PID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosKin{"Kinematics", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histos4piKin{"Four-Pion-Kinematics", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosMCtruth{"MC-Truth", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Debugging
  Configurable<bool> debugMode{"debugMode", false, "Enable Debug Mode"};
  // Configurable Event parameters
  Configurable<int> ifUPC{"ifUPC", 1, "Enable UPC reconstruction only"};
  Configurable<float> vZCut{"vZCut", 10., "Vertex Cut"};
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 50., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> zdcCut{"zdcCut", 1e6, "ZDC threshold"};
  Configurable<float> zdcMaxAmp{"zdcMaxAmp", 0, "ZDC max amplitude to be 0n"};
  Configurable<float> zdcMaxTime{"zdcMaxTime", 2, "ZDC max time in ns"};
  Configurable<std::string> neutronClass{"neutronClass", "XnXn", "Neutron class for ZDCs"};
  Configurable<uint16_t> numPVContrib{"numPVContrib", 4, "Number of PV Contributors"};
  Configurable<int> sbpCut{"sbpCut", 1, "Sbp"};
  Configurable<int> itsROFbCut{"itsROFbCut", 1, "itsROFbCut"};
  Configurable<int> vtxITSTPCcut{"vtxITSTPCcut", 1, "vtxITSTPCcut"};
  Configurable<int> tfbCut{"tfbCut", 1, "tfbCut"};
  // Configurable Track parameters
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use Only PV tracks"};
  Configurable<float> pTcut{"pTcut", 0.15, "Track Pt"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0, "dcaXY cut"};
  Configurable<float> dcaZcut{"dcaZcut", 2, "dcaZ cut"};
  Configurable<bool> useITStracksOnly{"useITStracksOnly", true, "only use tracks with hit in ITS"};
  Configurable<bool> useTPCtracksOnly{"useTPCtracksOnly", true, "only use tracks with hit in TPC"};
  Configurable<float> itsChi2NClsCut{"itsChi2NClsCut", 36, "ITS Chi2NCls"};
  Configurable<float> tpcChi2NClsCut{"tpcChi2NClsCut", 4.0, "TPC Chi2NCls"};
  Configurable<int> tpcNClsCrossedRowsCut{"tpcNClsCrossedRowsCut", 70, "Min TPC Findable Clusters"};
  // Configurable PID parameters
  Configurable<bool> ifCircularNSigmaCut{"ifCircularNSigmaCut", true, "Use circular nsigma cut for PID"};
  Configurable<bool> useTOF{"useTOF", true, "if track has TOF use TOF"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 5, "TPC cut"};
  Configurable<float> nSigmaTOFcut{"nSigmaTOFcut", 5, "TOF cut"};
  // Axis Configurations
  ConfigurableAxis pTAxis{"pTAxis", {1000, 0, 1}, "Axis for pT histograms"};
  ConfigurableAxis etaAxis{"etaAxis", {1000, -1.1, 1.1}, "Axis for Eta histograms"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {1000, -2.5, 2.5}, "Axis for Rapidity histograms"};
  ConfigurableAxis invMassAxis{"invMassAxis", {1000, 1, 2.5}, "Axis for Phi histograms"};
  ConfigurableAxis phiAxis{"phiAxis", {360, -1 * o2::constants::math::PI, o2::constants::math::PI}, "Axis for Phi histograms"};
  ConfigurableAxis cosThetaAxis{"cosThetaAxis", {360, -1, 1}, "Axis for cos Theta histograms"};

  void init(InitContext const&)
  {
    // QA plots: Event and Track Counter
    histosDataCounter.add("EventsCounts_vs_runNo", "Event Counter Run by Run; Run Number; Number of Events", kTH2F, {{113, 0, 113}, {14, 0, 14}});
    histosDataCounter.add("TracksCounts_vs_runNo", "Track Counter Run by Run; Run Number; Number of Track", kTH2F, {{113, 0, 113}, {13, 0, 13}});
    // QA plots: event selection-selected events
    histosQA.add("Events/selected/UPCmode", "UPC mode; Events", kTH1F, {{5, 0, 5}});
    histosQA.add("Events/selected/GapSide", "Gap Side;Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/TrueGapSide", "True Gap Side; True Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/isCBTOk", "isCBTOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/isCBTHadronOk", "isCBTHadronOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/isCBTZdcOk", "isCBTZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/isCBTHadronZdcOk", "isCBTHadronZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/selected/FT0A", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosQA.add("Events/selected/FT0C", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosQA.add("Events/selected/FV0A", "V0A amplitude", kTH1F, {{100, 0.0, 100}});
    histosQA.add("Events/selected/ZDC", "; ZDC A; ZDC C; time ZDC A [ns]; time ZDC C [ns]", kTHnSparseF, {{200, -10, 1000}, {200, -10, 1000}, {400, -10, 50}, {400, -10, 10}});
    histosQA.add("Events/selected/FDDA", "FDD A signal; FDD A signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosQA.add("Events/selected/FDDC", "FDD C signal; FDD C signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosQA.add("Events/selected/vertexX", "Vertex X; Vertex X [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosQA.add("Events/selected/vertexY", "Vertex Y; Vertex Y [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosQA.add("Events/selected/vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{2000, -15, 15}});
    histosQA.add("Events/selected/occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{20000, 0, 20000}});
    // QA plots: event selection-4 pion events
    histosQA.add("Events/4pion/UPCmode", "UPC mode; Events", kTH1F, {{5, 0, 5}});
    histosQA.add("Events/4pion/GapSide", "Gap Side;Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/TrueGapSide", "True Gap Side; True Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/isCBTOk", "isCBTOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/isCBTHadronOk", "isCBTHadronOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/isCBTZdcOk", "isCBTZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/isCBTHadronZdcOk", "isCBTHadronZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosQA.add("Events/4pion/FT0A", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosQA.add("Events/4pion/FT0C", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosQA.add("Events/4pion/FV0A", "V0A amplitude", kTH1F, {{100, 0.0, 100}});
    histosQA.add("Events/4pion/ZDC", "; ZDC A; ZDC C; time ZDC A; time ZDC C", kTHnSparseF, {{200, -10, 1000}, {200, -10, 1000}, {400, -10, 50}, {400, -10, 10}});
    histosQA.add("Events/4pion/FDDA", "FDD A signal; FDD A signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosQA.add("Events/4pion/FDDC", "FDD C signal; FDD C signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosQA.add("Events/4pion/vertexX", "Vertex X; Vertex X [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosQA.add("Events/4pion/vertexY", "Vertex Y; Vertex Y [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosQA.add("Events/4pion/vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{2000, -15, 15}});
    histosQA.add("Events/4pion/occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{20000, 0, 20000}});
    // QA plots: All tracks in selected events
    histosQA.add("Tracks/all/isPVcontributor", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{3, 0, 3}});
    histosQA.add("Tracks/all/dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/all/dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/all/itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/all/itsChi2", "ITS Chi2; ITS Chi2; Counts", kTH1F, {{500, 0, 50}});
    histosQA.add("Tracks/all/tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 10}});
    histosQA.add("Tracks/all/tpcNClsCrossedRows", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: Selected tracks in selected events
    histosQA.add("Tracks/selected/isPVcontributor", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{3, 0, 3}});
    histosQA.add("Tracks/selected/dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/selected/dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/selected/itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/selected/itsChi2", "ITS Chi2; ITS Chi2; Counts", kTH1F, {{500, 0, 50}});
    histosQA.add("Tracks/selected/tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/selected/tpcNClsCrossedRows", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: Pion tracks in selected events
    histosQA.add("Tracks/pions/isPVcontributor", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{3, 0, 3}});
    histosQA.add("Tracks/pions/dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/pions/dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/pions/itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/pions/itsChi2", "ITS Chi2; ITS Chi2; Counts", kTH1F, {{500, 0, 50}});
    histosQA.add("Tracks/pions/tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/pions/tpcNClsCrossedRows", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: Pion tracks from 4pi in selected events
    histosQA.add("Tracks/pions-from-4pi/isPVcontributor", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{3, 0, 3}});
    histosQA.add("Tracks/pions-from-4pi/dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/pions-from-4pi/dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosQA.add("Tracks/pions-from-4pi/itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/pions-from-4pi/itsChi2", "ITS Chi2; ITS Chi2; Counts", kTH1F, {{500, 0, 50}});
    histosQA.add("Tracks/pions-from-4pi/tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosQA.add("Tracks/pions-from-4pi/tpcNClsCrossedRows", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: PID- All tracks
    histosPID.add("all/tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 3}, {5000, 0.0, 600}});
    histosPID.add("all/tpcNSigmaPi", "TPC nSigma Pion for all tracks in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tpcNSigmaKa", "TPC nSigma Kaon for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tpcNSigmaPr", "TPC nSigma Proton for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tpcNSigmaEl", "TPC nSigma Electron for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tpcNSigmaMu", "TPC nSigma Muon for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tofBeta", "TOF beta vs p ; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {1500, 0.0, 1.5}});
    histosPID.add("all/tofNSigmaPi", "TOF nSigma Pion for all tracks in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tofNSigmaKa", "TOF nSigma Kaon for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tofNSigmaPr", "TOF nSigma Proton for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tofNSigmaEl", "TOF nSigma Electron for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("all/tofNSigmaMu", "TOF nSigma Muon for all tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    // QA plots: PID- Selected tracks
    histosPID.add("selected/tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 3}, {5000, 0.0, 600.0}});
    histosPID.add("selected/tpcNSigmaPi", "TPC nSigma Pion for all selected tracks in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tpcNSigmaKa", "TPC nSigma Kaon for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tpcNSigmaPr", "TPC nSigma Proton for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tpcNSigmaEl", "TPC nSigma Electron for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tpcNSigmaMu", "TPC nSigma Muon for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 2}, {1500, 0.0, 1.5}});
    histosPID.add("selected/tofNSigmaPi", "TOF nSigma Pion for all selected tracks in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tofNSigmaKa", "TOF nSigma Kaon for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tofNSigmaPr", "TOF nSigma Proton for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tofNSigmaEl", "TOF nSigma Electron for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("selected/tofNSigmaMu", "TOF nSigma Muon for all selected tracks in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    // QA plots: PID- Pion tracks
    histosPID.add("pions/tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 3}, {5000, 0.0, 600.0}});
    histosPID.add("pions/tpcNSigmaPi", "TPC nSigma Pion for all selected pions in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tpcNSigmaKa", "TPC nSigma Kaon for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tpcNSigmaPr", "TPC nSigma Proton for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tpcNSigmaEl", "TPC nSigma Electron for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tpcNSigmaMu", "TPC nSigma Muon for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {1500, 0.0, 1.5}});
    histosPID.add("pions/tofNSigmaPi", "TOF nSigma Pion for all selected pions in selected events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tofNSigmaKa", "TOF nSigma Kaon for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tofNSigmaPr", "TOF nSigma Proton for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tofNSigmaEl", "TOF nSigma Electron for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions/tofNSigmaMu", "TOF nSigma Muon for all selected pions in selected events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    // QA plots: PID- Pion tracks from 4pi events
    histosPID.add("pions-from-4pi/tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 3}, {5000, 0.0, 600.0}});
    histosPID.add("pions-from-4pi/tpcNSigmaPi", "TPC nSigma Pion for all pions from 4-pi events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tpcNSigmaKa", "TPC nSigma Kaon for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tpcNSigmaPr", "TPC nSigma Proton for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tpcNSigmaEl", "TPC nSigma Electron for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tpcNSigmaMu", "TPC nSigma Muon for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {1500, 0.0, 1.5}});
    histosPID.add("pions-from-4pi/tofNSigmaPi", "TOF nSigma Pion for all pions from 4-pi events; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tofNSigmaKa", "TOF nSigma Kaon for all pions from 4-pi eventsn; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tofNSigmaPr", "TOF nSigma Proton for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tofNSigmaEl", "TOF nSigma Electron for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosPID.add("pions-from-4pi/tofNSigmaMu", "TOF nSigma for all pions from 4-pi events; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    // Kinematics for all particles
    histosKin.add("all", ";pT [GeV/c]; #eta;#varphi", kTH3F, {pTAxis, etaAxis, phiAxis});
    histosKin.add("selected", ";pT [GeV/c]; #eta;#varphi", kTH3F, {pTAxis, etaAxis, phiAxis});
    histosKin.add("pions", ";pT [GeV/c]; #eta;#varphi", kTH3F, {pTAxis, etaAxis, phiAxis});
    histosKin.add("pions-from-4pion", ";pT [GeV/c]; #eta;#varphi;y ", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis});
    // Rho Prime Kinematics
    histos4piKin.add("two-pion", ";p_{T}^{4#pi} [GeV/c] ;m_{#pi^{+}#pi^{-}} [GeV/c^2];m_{#pi^{+}#pi^{-}} [GeV/c^2];m_{#pi^{+}#pi^{-}} [GeV/c^2];m_{#pi^{+}#pi^{-}} [GeV/c^2];m_{4#pi} [GeV/c^{2}]", kTHnSparseF, {{100, 0, 2}, {100, 0, 2}, {100, 0, 2}, {100, 0, 2}, invMassAxis});
    histos4piKin.add("zero-charge", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}]; Collin-Soper cos(#theta); Collin-Soper #varphi [rad];Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, cosThetaAxis, phiAxis, {113, 0, 113}});
    histos4piKin.add("non-zero-charge", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}];Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    histos4piKin.add("3piMinus-1piPlus", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}]; Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    histos4piKin.add("3piPlus-1piMinus", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}]; Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    histos4piKin.add("4piPlus", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}]; Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    histos4piKin.add("4piMinus", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}]; Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    // MC truth
    histosMCtruth.add("4-pi-pions", ";pT [GeV/c]; #eta;#varphi;y ", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, {113, 0, 113}});
    histosMCtruth.add("Four-pion", ";pT [GeV/c]; #eta; #varphi [rad];y; m_{4#pi} [GeV/c^{2}];Run Number", kTHnSparseF, {pTAxis, etaAxis, phiAxis, rapidityAxis, invMassAxis, {113, 0, 113}});
    //_______________________________________________________________________________________________________________________________________________
    setHistBinLabels();
    if (debugMode) {
      histosDataCounter.print();
      histosQA.print();
      histosPID.print();
      histosKin.print();
      histos4piKin.print();
      histosMCtruth.print();
    }
  } // End of init function

  //---------------------------------------------------------------------------------------------------------------------------------------------
  // Event Cuts
  Filter vertexZcut = (nabs(o2::aod::collision::posZ) <= vZCut);
  Filter numPVcontributorsCut = (o2::aod::collision::numContrib == numPVContrib);
  Filter fitcuts = (o2::aod::udcollision::totalFV0AmplitudeA <= fv0Cut) && (o2::aod::udcollision::totalFT0AmplitudeA <= ft0aCut) && (o2::aod::udcollision::totalFT0AmplitudeC <= ft0cCut);
  Filter zdcCuts = (o2::aod::udzdc::energyCommonZNA <= zdcCut) && (o2::aod::udzdc::energyCommonZNC <= zdcCut);
  Filter bcSelectionCuts = (o2::aod::udcollision::sbp == sbpCut) && (o2::aod::udcollision::itsROFb == itsROFbCut) && (o2::aod::udcollision::vtxITSTPC == vtxITSTPCcut) && (o2::aod::udcollision::tfb == tfbCut);
  // Track Cuts
  Filter tpcchi2nclsFilter = o2::aod::track::tpcChi2NCl <= tpcChi2NClsCut;
  Filter itschi2nclsFilter = o2::aod::track::itsChi2NCl <= itsChi2NClsCut;
  //---------------------------------------------------------------------------------------------------------------------------------------------

  using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;

  void processData(soa::Filtered<UDCollisions>::iterator const& collision, soa::Filtered<UDtracks> const& tracks)
  {

    // Check if the Event is reconstructed in UPC mode and RCT flag
    if ((collision.flags() != ifUPC) || (!sgSelector.isCBTHadronOk(collision))) {
      return;
    }

    int runIndex = getRunNumberIndex(collision.runNumber());

    histosQA.fill(HIST("Events/selected/UPCmode"), collision.flags());
    histosQA.fill(HIST("Events/selected/GapSide"), collision.gapSide());
    histosQA.fill(HIST("Events/selected/TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
    histosQA.fill(HIST("Events/selected/isCBTOk"), sgSelector.isCBTOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTHadronOk"), sgSelector.isCBTHadronOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTZdcOk"), sgSelector.isCBTZdcOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTHadronZdcOk"), sgSelector.isCBTHadronZdcOk(collision));
    histosQA.fill(HIST("Events/selected/vertexX"), collision.posX());
    histosQA.fill(HIST("Events/selected/vertexY"), collision.posY());
    histosQA.fill(HIST("Events/selected/vertexZ"), collision.posZ());
    histosQA.fill(HIST("Events/selected/occupancy"), collision.occupancyInTime());
    histosQA.fill(HIST("Events/selected/FV0A"), collision.totalFV0AmplitudeA());
    histosQA.fill(HIST("Events/selected/FT0A"), collision.totalFT0AmplitudeA());
    histosQA.fill(HIST("Events/selected/FT0C"), collision.totalFT0AmplitudeC());
    histosQA.fill(HIST("Events/selected/ZDC"), collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC());
    histosQA.fill(HIST("Events/selected/FDDA"), collision.totalFDDAmplitudeA());
    histosQA.fill(HIST("Events/selected/FDDC"), collision.totalFDDAmplitudeC());

    std::vector<decltype(tracks.begin())> selectedPionTracks;
    std::vector<decltype(tracks.begin())> selectedPionPlusTracks;
    std::vector<decltype(tracks.begin())> selectedPionMinusTracks;

    for (const auto& t0 : tracks) {

      PxPyPzMVector tVector(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);

      // QA-Tracks before selection
      histosQA.fill(HIST("Tracks/all/isPVcontributor"), t0.isPVContributor());
      histosQA.fill(HIST("Tracks/all/dcaXY"), t0.dcaXY());
      histosQA.fill(HIST("Tracks/all/dcaZ"), t0.dcaZ());
      histosQA.fill(HIST("Tracks/all/itsChi2NCl"), t0.itsChi2NCl());
      histosQA.fill(HIST("Tracks/all/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
      histosQA.fill(HIST("Tracks/all/tpcChi2NCl"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/all/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

      // PID before track selection
      histosPID.fill(HIST("all/tpcSignal"), tVector.P(), t0.tpcSignal());
      histosPID.fill(HIST("all/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
      histosPID.fill(HIST("all/tofBeta"), tVector.P(), t0.beta());
      histosPID.fill(HIST("all/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

      // Kinematics for all particles before selection
      histosKin.fill(HIST("all"), tVector.Pt(), tVector.Eta(), tVector.Phi());

      // Selecting good tracks
      if (!isSelectedTrack(t0, pTcut, etaCut, dcaXYcut, dcaZcut, useITStracksOnly, useTPCtracksOnly, itsChi2NClsCut, tpcChi2NClsCut, tpcNClsCrossedRowsCut)) {
        continue;
      }

      // QA-Tracks after selection
      histosQA.fill(HIST("Tracks/selected/isPVcontributor"), t0.isPVContributor());
      histosQA.fill(HIST("Tracks/selected/dcaXY"), t0.dcaXY());
      histosQA.fill(HIST("Tracks/selected/dcaZ"), t0.dcaZ());
      histosQA.fill(HIST("Tracks/selected/itsChi2NCl"), t0.itsChi2NCl());
      histosQA.fill(HIST("Tracks/selected/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
      histosQA.fill(HIST("Tracks/selected/tpcChi2NCl"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/selected/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

      // PID after track selection before selecting pions
      histosPID.fill(HIST("selected/tpcSignal"), tVector.P(), t0.tpcSignal());
      histosPID.fill(HIST("selected/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
      histosPID.fill(HIST("selected/tofBeta"), tVector.P(), t0.beta());
      histosPID.fill(HIST("selected/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

      // Kinematics for all particles after track selection before selecting pions
      histosKin.fill(HIST("selected"), tVector.Pt(), tVector.Eta(), tVector.Phi());

      if (ifPion(t0, useTOF, nSigmaTPCcut, nSigmaTOFcut, ifCircularNSigmaCut)) {

        selectedPionTracks.push_back(t0);

        // QA-Tracks after selecting pions
        histosQA.fill(HIST("Tracks/pions/isPVcontributor"), t0.isPVContributor());
        histosQA.fill(HIST("Tracks/pions/dcaXY"), t0.dcaXY());
        histosQA.fill(HIST("Tracks/pions/dcaZ"), t0.dcaZ());
        histosQA.fill(HIST("Tracks/pions/itsChi2NCl"), t0.itsChi2NCl());
        histosQA.fill(HIST("Tracks/pions/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
        histosQA.fill(HIST("Tracks/pions/tpcChi2NCl"), t0.tpcChi2NCl());
        histosQA.fill(HIST("Tracks/pions/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

        // PID after selecting pions
        histosPID.fill(HIST("pions/tpcSignal"), tVector.P(), t0.tpcSignal());
        histosPID.fill(HIST("pions/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
        histosPID.fill(HIST("pions/tofBeta"), tVector.P(), t0.beta());
        histosPID.fill(HIST("pions/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

        // Kinematics for pions
        histosKin.fill(HIST("pions"), tVector.Pt(), tVector.Eta(), tVector.Phi());

        if (t0.sign() == 1) {
          selectedPionPlusTracks.push_back(t0);
        }
        if (t0.sign() == -1) {
          selectedPionMinusTracks.push_back(t0);
        }
      } // End of Selection PID Pion
    } // End of loop over tracks

    int numSelectedPionTracks = static_cast<int>(selectedPionTracks.size());
    int numPiPlusTracks = static_cast<int>(selectedPionPlusTracks.size());
    int numPionMinusTracks = static_cast<int>(selectedPionMinusTracks.size());

    // event should have exactly 4 pions
    if (numSelectedPionTracks != four) {
      return;
    }

    // Selecting Events with net charge = 0
    if (numPionMinusTracks == two && numPiPlusTracks == two) {

      // QA-Events-4pion
      histosQA.fill(HIST("Events/4pion/UPCmode"), collision.flags());
      histosQA.fill(HIST("Events/4pion/GapSide"), collision.gapSide());
      histosQA.fill(HIST("Events/4pion/TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
      histosQA.fill(HIST("Events/4pion/isCBTOk"), sgSelector.isCBTOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTHadronOk"), sgSelector.isCBTHadronOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTZdcOk"), sgSelector.isCBTZdcOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTHadronZdcOk"), sgSelector.isCBTHadronZdcOk(collision));
      histosQA.fill(HIST("Events/4pion/vertexX"), collision.posX());
      histosQA.fill(HIST("Events/4pion/vertexY"), collision.posY());
      histosQA.fill(HIST("Events/4pion/vertexZ"), collision.posZ());
      histosQA.fill(HIST("Events/4pion/occupancy"), collision.occupancyInTime());
      histosQA.fill(HIST("Events/4pion/FV0A"), collision.totalFV0AmplitudeA());
      histosQA.fill(HIST("Events/4pion/FT0A"), collision.totalFT0AmplitudeA());
      histosQA.fill(HIST("Events/4pion/FT0C"), collision.totalFT0AmplitudeC());
      histosQA.fill(HIST("Events/4pion/ZDC"), collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC());
      histosQA.fill(HIST("Events/4pion/FDDA"), collision.totalFDDAmplitudeA());
      histosQA.fill(HIST("Events/4pion/FDDC"), collision.totalFDDAmplitudeC());

      for (int i = 0; i < four; i++) {
        PxPyPzMVector tVector(selectedPionTracks[i].px(), selectedPionTracks[i].py(), selectedPionTracks[i].pz(), o2::constants::physics::MassPionCharged);
        // Tracks QA for all four pions
        histosQA.fill(HIST("Tracks/pions-from-4pi/isPVcontributor"), selectedPionTracks[i].isPVContributor());
        histosQA.fill(HIST("Tracks/pions-from-4pi/dcaXY"), selectedPionTracks[i].dcaXY());
        histosQA.fill(HIST("Tracks/pions-from-4pi/dcaZ"), selectedPionTracks[i].dcaZ());
        histosQA.fill(HIST("Tracks/pions-from-4pi/itsChi2NCl"), selectedPionTracks[i].itsChi2NCl());
        histosQA.fill(HIST("Tracks/pions-from-4pi/itsChi2"), selectedPionTracks[i].itsChi2NCl() * selectedPionTracks[i].itsNCls());
        histosQA.fill(HIST("Tracks/pions-from-4pi/tpcChi2NCl"), selectedPionTracks[i].tpcChi2NCl());
        histosQA.fill(HIST("Tracks/pions-from-4pi/tpcNClsCrossedRows"), selectedPionTracks[i].tpcNClsCrossedRows());
        // PID for all four pions
        histosPID.fill(HIST("pions-from-4pi/tpcSignal"), tVector.P(), selectedPionTracks[i].tpcSignal());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaPi"), selectedPionTracks[i].tpcNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaKa"), selectedPionTracks[i].tpcNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaPr"), selectedPionTracks[i].tpcNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaEl"), selectedPionTracks[i].tpcNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaMu"), selectedPionTracks[i].tpcNSigmaMu(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofBeta"), tVector.P(), selectedPionTracks[i].beta());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaPi"), selectedPionTracks[i].tofNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaKa"), selectedPionTracks[i].tofNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaPr"), selectedPionTracks[i].tofNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaEl"), selectedPionTracks[i].tofNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaMu"), selectedPionTracks[i].tofNSigmaMu(), tVector.Pt());
      }

      PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

      // Kinematics for pions from 4 pion events
      histosKin.fill(HIST("pions-from-4pion"), p1.Pt(), p1.Eta(), p1.Phi(), p1.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p2.Pt(), p2.Eta(), p2.Phi(), p2.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p3.Pt(), p3.Eta(), p3.Phi(), p3.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p4.Pt(), p4.Eta(), p4.Phi(), p4.Rapidity());

      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
      PxPyPzMVector p13 = p1 + p3;
      PxPyPzMVector p14 = p1 + p4;
      PxPyPzMVector p23 = p2 + p3;
      PxPyPzMVector p24 = p2 + p4;

      // Two Pion Mass combinations
      histos4piKin.fill(HIST("two-pion"), p1234.Pt(), p13.M(), p14.M(), p23.M(), p24.M(), p1234.M());

      double fourPiPhiPair1 = collinSoperPhi(p13, p1234);
      double fourPiPhiPair2 = collinSoperPhi(p14, p1234);
      double fourPiPhiPair3 = collinSoperPhi(p23, p1234);
      double fourPiPhiPair4 = collinSoperPhi(p24, p1234);

      double fourPiCosThetaPair1 = collinSoperCosTheta(p13, p1234);
      double fourPiCosThetaPair2 = collinSoperCosTheta(p14, p1234);
      double fourPiCosThetaPair3 = collinSoperCosTheta(p23, p1234);
      double fourPiCosThetaPair4 = collinSoperCosTheta(p24, p1234);

      double mDiff13 = std::abs((p13.M() - mRho0));
      double mDiff14 = std::abs((p14.M() - mRho0));
      double mDiff23 = std::abs((p23.M() - mRho0));
      double mDiff24 = std::abs((p24.M() - mRho0));
      if ((mDiff13 < mDiff14) && (mDiff13 < mDiff23) && (mDiff13 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair1, fourPiPhiPair1, runIndex);
      } else if ((mDiff14 < mDiff13) && (mDiff14 < mDiff23) && (mDiff14 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair2, fourPiPhiPair2, runIndex);
      } else if ((mDiff23 < mDiff13) && (mDiff23 < mDiff14) && (mDiff23 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair3, fourPiPhiPair3, runIndex);
      } else if ((mDiff24 < mDiff13) && (mDiff24 < mDiff14) && (mDiff24 < mDiff23)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair4, fourPiPhiPair4, runIndex);
      }
    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTracks != two && numPiPlusTracks != two) {
      PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
      // Kinematics for 4 pion system from non 0 charge events
      if (numPionMinusTracks == three && numPiPlusTracks == one) {
        histos4piKin.fill(HIST("3piMinus-1piPlus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == one && numPiPlusTracks == three) {
        histos4piKin.fill(HIST("3piPlus-1piMinus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == four && numPiPlusTracks == zero) {
        histos4piKin.fill(HIST("4piMinus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == zero && numPiPlusTracks == four) {
        histos4piKin.fill(HIST("4piPlus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      }
      histos4piKin.fill(HIST("non-zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
    } // End of Analysis for non 0 charge events
  } // End of 4 Pion Analysis Process function for Pass5 Data

  void processEventCounter(UDCollisions::iterator const& collision)
  {
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 0);
    // RCT flag
    if (!sgSelector.isCBTHadronZdcOk(collision)) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 1);
    // UPC mode
    if (collision.flags() != ifUPC) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 2);
    // vtxITSTPC
    if (collision.vtxITSTPC() != vtxITSTPCcut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 3);
    // sbp
    if (collision.sbp() != sbpCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 4);
    // itsROFb
    if (collision.itsROFb() != itsROFbCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 5);
    // tfb
    if (collision.tfb() != tfbCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 6);
    // FT0A
    if (collision.totalFT0AmplitudeA() > ft0aCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 7);
    // FT0C
    if (collision.totalFT0AmplitudeC() > ft0cCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 8);
    // FV0A
    if (collision.totalFV0AmplitudeA() > fv0Cut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 9);
    // ZDC
    if (!neutronClassSelection(collision)) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 10);
    // numContributors
    if (collision.numContrib() != numPVContrib) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 11);
    // vertexZ
    if (std::abs(collision.posZ()) > vZCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 12);
  } // End of processCounter function

  void processTrackCounter(soa::Filtered<UDCollisions>::iterator const& collision, UDtracks const& tracks)
  {

    int runIndex = getRunNumberIndex(collision.runNumber());

    // Check if the Event is reconstructed in UPC mode
    if ((collision.flags() != ifUPC) || (!sgSelector.isCBTHadronZdcOk(collision))) {
      return;
    }

    for (const auto& track : tracks) {

      // total tracks
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 0);

      PxPyPzMVector trackVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);

      // pt cut
      if (trackVector.Pt() < pTcut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 1);

      // eta cut
      if (std::abs(trackVector.Eta()) > etaCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 2);

      // DCA Z cut
      if (std::abs(track.dcaZ()) > dcaZcut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 3);

      // DCA XY cut
      float maxDCAxy = 0.0105 + 0.035 / std::pow(trackVector.Pt(), 1.1);
      if (dcaXYcut == 0 && (std::fabs(track.dcaXY()) > maxDCAxy)) {
        continue;
      } else if (dcaXYcut != 0 && (std::fabs(track.dcaXY()) > dcaXYcut)) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 4);

      // ITS Track only
      if (useITStracksOnly && !track.hasITS()) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 5);

      // TPC Track only
      if (useTPCtracksOnly && !track.hasTPC()) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 6);

      // ITS Chi2 N Clusters cut
      if (track.hasITS() && track.itsChi2NCl() > itsChi2NClsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 7);

      // TPC Chi2 N Clusters cut
      if (track.hasTPC() && track.tpcChi2NCl() > tpcChi2NClsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 8);

      // TPC N Clusters Findable cut
      if (track.hasTPC() && track.tpcNClsCrossedRows() < tpcNClsCrossedRowsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 9);

      // Selection PID Pion
      if (ifPion(track, useTOF, nSigmaTPCcut, nSigmaTOFcut, ifCircularNSigmaCut)) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 10);

      // is PV contributor
      if (track.isPVContributor() != useOnlyPVtracks) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 11);

    } // End of loop over tracks
  } // End of processCounter function

  using MCtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;
  using MCCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDMcCollsLabels>;

  void processMCgen(aod::UDMcCollisions::iterator const&, aod::UDMcParticles const& mcParticles, aod::BCs const& bcs)
  {

    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    int runIndex = getRunNumberIndex(bc.runNumber());

    for (const auto& particle : mcParticles) {
      PxPyPzMVector p1234;
      if ((particle.pdgCode() != rhoPrime) || (particle.daughters_as<aod::UDMcParticles>().size() != four)) {
        continue;
      }
      for (const auto& daughter : particle.daughters_as<aod::UDMcParticles>()) {
        PxPyPzMVector dVector(daughter.px(), daughter.py(), daughter.pz(), o2::constants::physics::MassPionCharged);
        if (daughter.pdgCode() == PDG_t::kPiPlus) {
          histosMCtruth.fill(HIST("4-pi-pions"), dVector.Pt(), dVector.Eta(), dVector.Phi(), dVector.Rapidity(), runIndex);
          p1234 = p1234 + dVector;
        }
        if (daughter.pdgCode() == PDG_t::kPiMinus) {
          histosMCtruth.fill(HIST("4-pi-pions"), dVector.Pt(), dVector.Eta(), dVector.Phi(), dVector.Rapidity(), runIndex);
          p1234 = p1234 + dVector;
        }
      } // End of loop over daughters
      histosMCtruth.fill(HIST("Four-pion"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
    } // End of loop over MC particles
  } // End of processMCgen function

  void processMCrec(soa::Filtered<MCCollisions>::iterator const& collision, soa::Filtered<MCtracks> const& tracks)
  {

    // Check if the Event is reconstructed in UPC mode and RCT flag
    if ((collision.flags() != ifUPC) || (!sgSelector.isCBTHadronOk(collision)) || (!collision.has_udMcCollision())) {
      return;
    }

    int runIndex = getRunNumberIndex(collision.runNumber());

    histosQA.fill(HIST("Events/selected/UPCmode"), collision.flags());
    histosQA.fill(HIST("Events/selected/GapSide"), collision.gapSide());
    histosQA.fill(HIST("Events/selected/TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
    histosQA.fill(HIST("Events/selected/isCBTOk"), sgSelector.isCBTOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTHadronOk"), sgSelector.isCBTHadronOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTZdcOk"), sgSelector.isCBTZdcOk(collision));
    histosQA.fill(HIST("Events/selected/isCBTHadronZdcOk"), sgSelector.isCBTHadronZdcOk(collision));
    histosQA.fill(HIST("Events/selected/vertexX"), collision.posX());
    histosQA.fill(HIST("Events/selected/vertexY"), collision.posY());
    histosQA.fill(HIST("Events/selected/vertexZ"), collision.posZ());
    histosQA.fill(HIST("Events/selected/occupancy"), collision.occupancyInTime());
    histosQA.fill(HIST("Events/selected/FV0A"), collision.totalFV0AmplitudeA());
    histosQA.fill(HIST("Events/selected/FT0A"), collision.totalFT0AmplitudeA());
    histosQA.fill(HIST("Events/selected/FT0C"), collision.totalFT0AmplitudeC());
    histosQA.fill(HIST("Events/selected/ZDC"), collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC());
    histosQA.fill(HIST("Events/selected/FDDA"), collision.totalFDDAmplitudeA());
    histosQA.fill(HIST("Events/selected/FDDC"), collision.totalFDDAmplitudeC());

    std::vector<decltype(tracks.begin())> selectedPionTracks;
    std::vector<decltype(tracks.begin())> selectedPionPlusTracks;
    std::vector<decltype(tracks.begin())> selectedPionMinusTracks;

    for (const auto& t0 : tracks) {

      PxPyPzMVector tVector(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);

      // QA-Tracks before selection
      histosQA.fill(HIST("Tracks/all/isPVcontributor"), t0.isPVContributor());
      histosQA.fill(HIST("Tracks/all/dcaXY"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/all/dcaZ"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/all/itsChi2NCl"), t0.itsChi2NCl());
      histosQA.fill(HIST("Tracks/all/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
      histosQA.fill(HIST("Tracks/all/tpcChi2NCl"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/all/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

      // PID before track selection
      histosPID.fill(HIST("all/tpcSignal"), tVector.P(), t0.tpcSignal());
      histosPID.fill(HIST("all/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("all/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
      histosPID.fill(HIST("all/tofBeta"), tVector.P(), t0.beta());
      histosPID.fill(HIST("all/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("all/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

      // Kinematics for all particles before selection
      histosKin.fill(HIST("all"), tVector.Pt(), tVector.Eta(), tVector.Phi());

      // Selecting good tracks
      if (!isSelectedTrack(t0, pTcut, etaCut, dcaXYcut, dcaZcut, useITStracksOnly, useTPCtracksOnly, itsChi2NClsCut, tpcChi2NClsCut, tpcNClsCrossedRowsCut)) {
        continue;
      }
      if (!t0.has_udMcParticle()) {
        continue;
      }

      // QA-Tracks after selection
      histosQA.fill(HIST("Tracks/selected/isPVcontributor"), t0.isPVContributor());
      histosQA.fill(HIST("Tracks/selected/dcaXY"), t0.dcaXY());
      histosQA.fill(HIST("Tracks/selected/dcaZ"), t0.dcaZ());
      histosQA.fill(HIST("Tracks/selected/itsChi2NCl"), t0.itsChi2NCl());
      histosQA.fill(HIST("Tracks/selected/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
      histosQA.fill(HIST("Tracks/selected/tpcChi2NCl"), t0.tpcChi2NCl());
      histosQA.fill(HIST("Tracks/selected/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

      // PID after track selection before selecting pions
      histosPID.fill(HIST("selected/tpcSignal"), tVector.P(), t0.tpcSignal());
      histosPID.fill(HIST("selected/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("selected/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
      histosPID.fill(HIST("selected/tofBeta"), tVector.P(), t0.beta());
      histosPID.fill(HIST("selected/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
      histosPID.fill(HIST("selected/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

      // Kinematics for all particles after track selection before selecting pions
      histosKin.fill(HIST("selected"), tVector.Pt(), tVector.Eta(), tVector.Phi());

      if (ifPion(t0, useTOF, nSigmaTPCcut, nSigmaTOFcut, ifCircularNSigmaCut)) {

        selectedPionTracks.push_back(t0);

        // QA-Tracks after selecting pions
        histosQA.fill(HIST("Tracks/pions/isPVcontributor"), t0.isPVContributor());
        histosQA.fill(HIST("Tracks/pions/dcaXY"), t0.dcaXY());
        histosQA.fill(HIST("Tracks/pions/dcaZ"), t0.dcaZ());
        histosQA.fill(HIST("Tracks/pions/itsChi2NCl"), t0.itsChi2NCl());
        histosQA.fill(HIST("Tracks/pions/itsChi2"), t0.itsChi2NCl() * t0.itsNCls());
        histosQA.fill(HIST("Tracks/pions/tpcChi2NCl"), t0.tpcChi2NCl());
        histosQA.fill(HIST("Tracks/pions/tpcNClsCrossedRows"), t0.tpcNClsCrossedRows());

        // PID after selecting pions
        histosPID.fill(HIST("pions/tpcSignal"), tVector.P(), t0.tpcSignal());
        histosPID.fill(HIST("pions/tpcNSigmaPi"), t0.tpcNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaKa"), t0.tpcNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaPr"), t0.tpcNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaEl"), t0.tpcNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions/tpcNSigmaMu"), t0.tpcNSigmaMu(), tVector.Pt());
        histosPID.fill(HIST("pions/tofBeta"), tVector.P(), t0.beta());
        histosPID.fill(HIST("pions/tofNSigmaPi"), t0.tofNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaKa"), t0.tofNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaPr"), t0.tofNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaEl"), t0.tofNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions/tofNSigmaMu"), t0.tofNSigmaMu(), tVector.Pt());

        // Kinematics for pions
        histosKin.fill(HIST("pions"), tVector.Pt(), tVector.Eta(), tVector.Phi());

        if (t0.sign() == 1) {
          selectedPionPlusTracks.push_back(t0);
        }
        if (t0.sign() == -1) {
          selectedPionMinusTracks.push_back(t0);
        }
      } // End of Selection PID Pion
    } // End of loop over tracks

    int numSelectedPionTracks = static_cast<int>(selectedPionTracks.size());
    int numPiPlusTracks = static_cast<int>(selectedPionPlusTracks.size());
    int numPionMinusTracks = static_cast<int>(selectedPionMinusTracks.size());

    // event should have exactly 4 pions
    if (numSelectedPionTracks != four) {
      return;
    }

    // Selecting Events with net charge = 0
    if (numPionMinusTracks == two && numPiPlusTracks == two) {

      // QA-Events-4pion
      histosQA.fill(HIST("Events/4pion/UPCmode"), collision.flags());
      histosQA.fill(HIST("Events/4pion/GapSide"), collision.gapSide());
      histosQA.fill(HIST("Events/4pion/TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
      histosQA.fill(HIST("Events/4pion/isCBTOk"), sgSelector.isCBTOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTHadronOk"), sgSelector.isCBTHadronOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTZdcOk"), sgSelector.isCBTZdcOk(collision));
      histosQA.fill(HIST("Events/4pion/isCBTHadronZdcOk"), sgSelector.isCBTHadronZdcOk(collision));
      histosQA.fill(HIST("Events/4pion/vertexX"), collision.posX());
      histosQA.fill(HIST("Events/4pion/vertexY"), collision.posY());
      histosQA.fill(HIST("Events/4pion/vertexZ"), collision.posZ());
      histosQA.fill(HIST("Events/4pion/occupancy"), collision.occupancyInTime());
      histosQA.fill(HIST("Events/4pion/FV0A"), collision.totalFV0AmplitudeA());
      histosQA.fill(HIST("Events/4pion/FT0A"), collision.totalFT0AmplitudeA());
      histosQA.fill(HIST("Events/4pion/FT0C"), collision.totalFT0AmplitudeC());
      histosQA.fill(HIST("Events/4pion/ZDC"), collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC());
      histosQA.fill(HIST("Events/4pion/FDDA"), collision.totalFDDAmplitudeA());
      histosQA.fill(HIST("Events/4pion/FDDC"), collision.totalFDDAmplitudeC());

      for (int i = 0; i < four; i++) {
        PxPyPzMVector tVector(selectedPionTracks[i].px(), selectedPionTracks[i].py(), selectedPionTracks[i].pz(), o2::constants::physics::MassPionCharged);
        // Tracks QA for all four pions
        histosQA.fill(HIST("Tracks/pions-from-4pi/isPVcontributor"), selectedPionTracks[i].isPVContributor());
        histosQA.fill(HIST("Tracks/pions-from-4pi/dcaXY"), selectedPionTracks[i].dcaXY());
        histosQA.fill(HIST("Tracks/pions-from-4pi/dcaZ"), selectedPionTracks[i].dcaZ());
        histosQA.fill(HIST("Tracks/pions-from-4pi/itsChi2NCl"), selectedPionTracks[i].itsChi2NCl());
        histosQA.fill(HIST("Tracks/pions-from-4pi/itsChi2"), selectedPionTracks[i].itsChi2NCl() * selectedPionTracks[i].itsNCls());
        histosQA.fill(HIST("Tracks/pions-from-4pi/tpcChi2NCl"), selectedPionTracks[i].tpcChi2NCl());
        histosQA.fill(HIST("Tracks/pions-from-4pi/tpcNClsCrossedRows"), selectedPionTracks[i].tpcNClsCrossedRows());
        // PID for all four pions
        histosPID.fill(HIST("pions-from-4pi/tpcSignal"), tVector.P(), selectedPionTracks[i].tpcSignal());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaPi"), selectedPionTracks[i].tpcNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaKa"), selectedPionTracks[i].tpcNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaPr"), selectedPionTracks[i].tpcNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaEl"), selectedPionTracks[i].tpcNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tpcNSigmaMu"), selectedPionTracks[i].tpcNSigmaMu(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofBeta"), tVector.P(), selectedPionTracks[i].beta());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaPi"), selectedPionTracks[i].tofNSigmaPi(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaKa"), selectedPionTracks[i].tofNSigmaKa(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaPr"), selectedPionTracks[i].tofNSigmaPr(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaEl"), selectedPionTracks[i].tofNSigmaEl(), tVector.Pt());
        histosPID.fill(HIST("pions-from-4pi/tofNSigmaMu"), selectedPionTracks[i].tofNSigmaMu(), tVector.Pt());
      }

      PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

      // Kinematics for pions from 4 pion events
      histosKin.fill(HIST("pions-from-4pion"), p1.Pt(), p1.Eta(), p1.Phi(), p1.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p2.Pt(), p2.Eta(), p2.Phi(), p2.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p3.Pt(), p3.Eta(), p3.Phi(), p3.Rapidity());
      histosKin.fill(HIST("pions-from-4pion"), p4.Pt(), p4.Eta(), p4.Phi(), p4.Rapidity());

      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
      PxPyPzMVector p13 = p1 + p3;
      PxPyPzMVector p14 = p1 + p4;
      PxPyPzMVector p23 = p2 + p3;
      PxPyPzMVector p24 = p2 + p4;

      // Two Pion Mass combinations
      histos4piKin.fill(HIST("two-pion"), p13.M(), p14.M(), p23.M(), p24.M(), p1234.M());

      double fourPiPhiPair1 = collinSoperPhi(p13, p1234);
      double fourPiPhiPair2 = collinSoperPhi(p14, p1234);
      double fourPiPhiPair3 = collinSoperPhi(p23, p1234);
      double fourPiPhiPair4 = collinSoperPhi(p24, p1234);

      double fourPiCosThetaPair1 = collinSoperCosTheta(p13, p1234);
      double fourPiCosThetaPair2 = collinSoperCosTheta(p14, p1234);
      double fourPiCosThetaPair3 = collinSoperCosTheta(p23, p1234);
      double fourPiCosThetaPair4 = collinSoperCosTheta(p24, p1234);

      double mDiff13 = std::abs((p13.M() - mRho0));
      double mDiff14 = std::abs((p14.M() - mRho0));
      double mDiff23 = std::abs((p23.M() - mRho0));
      double mDiff24 = std::abs((p24.M() - mRho0));

      if ((mDiff13 < mDiff14) && (mDiff13 < mDiff23) && (mDiff13 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair1, fourPiPhiPair1, runIndex);
      } else if ((mDiff14 < mDiff13) && (mDiff14 < mDiff23) && (mDiff14 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair2, fourPiPhiPair2, runIndex);
      } else if ((mDiff23 < mDiff13) && (mDiff23 < mDiff14) && (mDiff23 < mDiff24)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair3, fourPiPhiPair3, runIndex);
      } else if ((mDiff24 < mDiff13) && (mDiff24 < mDiff14) && (mDiff24 < mDiff23)) {
        histos4piKin.fill(HIST("zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), fourPiCosThetaPair4, fourPiPhiPair4, runIndex);
      }
    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTracks != two && numPiPlusTracks != two) {
      PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
      // Kinematics for 4 pion system from non 0 charge events
      if (numPionMinusTracks == three && numPiPlusTracks == one) {
        histos4piKin.fill(HIST("3piMinus-1piPlus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == one && numPiPlusTracks == three) {
        histos4piKin.fill(HIST("3piPlus-1piMinus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == four && numPiPlusTracks == zero) {
        histos4piKin.fill(HIST("4piMinus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      } else if (numPionMinusTracks == zero && numPiPlusTracks == four) {
        histos4piKin.fill(HIST("4piPlus"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
      }
      histos4piKin.fill(HIST("non-zero-charge"), p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(), runIndex);
    } // End of Analysis for non 0 charge events
  } // End of 4 Pion Analysis Process function for Pass5 MC

  void processEventCounterMC(MCCollisions::iterator const& collision)
  {

    // Check if the Event has MC labels
    if (!collision.has_udMcCollision()) {
      return;
    }

    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 0);
    // RCT flag
    if (!sgSelector.isCBTHadronZdcOk(collision)) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 1);
    // UPC mode
    if (collision.flags() != ifUPC) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 2);
    // vtxITSTPC
    if (collision.vtxITSTPC() != vtxITSTPCcut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 3);
    // sbp
    if (collision.sbp() != sbpCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 4);
    // itsROFb
    if (collision.itsROFb() != itsROFbCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 5);
    // tfb
    if (collision.tfb() != tfbCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 6);
    // FT0A
    if (collision.totalFT0AmplitudeA() > ft0aCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 7);
    // FT0C
    if (collision.totalFT0AmplitudeC() > ft0cCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 8);
    // FV0A
    if (collision.totalFV0AmplitudeA() > fv0Cut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 9);
    // ZDC
    if (collision.energyCommonZNA() > zdcCut || collision.energyCommonZNC() > zdcCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 10);
    // numContributors
    if (collision.numContrib() != numPVContrib) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 11);
    // vertexZ
    if (std::abs(collision.posZ()) > vZCut) {
      return;
    }
    histosDataCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 12);
  } // End of processCounter function

  void processTrackCounterMC(soa::Filtered<MCCollisions>::iterator const& collision, MCtracks const& tracks)
  {

    int runIndex = getRunNumberIndex(collision.runNumber());

    // Check if the Event is reconstructed in UPC mode
    if ((collision.flags() != ifUPC) || (!sgSelector.isCBTHadronZdcOk(collision))) {
      return;
    }

    for (const auto& track : tracks) {

      if (!track.has_udMcParticle()) {
        continue;
      }

      // total tracks
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 0);

      PxPyPzMVector trackVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);

      // pt cut
      if (trackVector.Pt() < pTcut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 1);

      // eta cut
      if (std::abs(trackVector.Eta()) > etaCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 2);

      // DCA Z cut
      if (std::abs(track.dcaZ()) > dcaZcut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 3);

      // DCA XY cut
      float maxDCAxy = 0.0105 + 0.035 / std::pow(trackVector.Pt(), 1.1);
      if (dcaXYcut == 0 && (std::fabs(track.dcaXY()) > maxDCAxy)) {
        continue;
      } else if (dcaXYcut != 0 && (std::fabs(track.dcaXY()) > dcaXYcut)) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 4);

      // ITS Track only
      if (useITStracksOnly && !track.hasITS()) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 5);

      // TPC Track only
      if (useTPCtracksOnly && !track.hasTPC()) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 6);

      // ITS Chi2 N Clusters cut
      if (track.hasITS() && track.itsChi2NCl() > itsChi2NClsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 7);

      // TPC Chi2 N Clusters cut
      if (track.hasTPC() && track.tpcChi2NCl() > tpcChi2NClsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 8);

      // TPC N Clusters Findable cut
      if (track.hasTPC() && track.tpcNClsCrossedRows() < tpcNClsCrossedRowsCut) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 9);

      // Selection PID Pion
      if (ifPion(track, useTOF, nSigmaTPCcut, nSigmaTOFcut, ifCircularNSigmaCut)) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 10);

      // is PV contributor
      if (track.isPVContributor() != useOnlyPVtracks) {
        continue;
      }
      histosDataCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 11);

    } // End of loop over tracks
  } // End of processCounter function

  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processData, "Data Analysis Function", true);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processEventCounter, "Event Counter Function", true);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processTrackCounter, "Track Counter Function", true);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processMCgen, "MC generated Analysis Function", false);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processMCrec, "MC reconstructed Analysis Function", false);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processEventCounterMC, "MC Event Counter Function", false);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processTrackCounterMC, "MC Track Counter Function", false);

  double collinSoperPhi(PxPyPzMVector twoPionVector, PxPyPzMVector fourPionVector)
  {
    PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target
    // Boost to center of mass frame
    Boost boosTo4PiCM{fourPionVector.BoostToCM()};
    XYZVectorF twoPionVectorCM{(boosTo4PiCM(twoPionVector).Vect()).Unit()};
    XYZVectorF beam1CM{(boosTo4PiCM(pProjCM).Vect()).Unit()};
    XYZVectorF beam2CM{(boosTo4PiCM(pTargCM).Vect()).Unit()};
    // Axes
    XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};
    XYZVectorF yaxisCS{(beam1CM.Cross(beam2CM)).Unit()};
    XYZVectorF xaxisCS{(yaxisCS.Cross(zaxisCS)).Unit()};
    double phi = std::atan2(yaxisCS.Dot(twoPionVectorCM), xaxisCS.Dot(twoPionVectorCM));
    return phi;
  }

  double collinSoperCosTheta(PxPyPzMVector twoPionVector, PxPyPzMVector fourPionVector)
  {
    PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target
    // Boost to center of mass frame
    Boost boosTo4PiCM{fourPionVector.BoostToCM()};
    XYZVectorF twoPionVectorCM{(boosTo4PiCM(twoPionVector).Vect()).Unit()};
    XYZVectorF beam1CM{(boosTo4PiCM(pProjCM).Vect()).Unit()};
    XYZVectorF beam2CM{(boosTo4PiCM(pTargCM).Vect()).Unit()};
    // Axes
    XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};
    double cosThetaCS = zaxisCS.Dot(twoPionVectorCM);
    return cosThetaCS;
  }

  template <typename C>
  bool neutronClassSelection(C const& coll)
  {

    bool aXn = coll.energyCommonZNA() > zdcMaxAmp && coll.timeZNA() < zdcMaxTime;
    bool a0n = coll.energyCommonZNA() <= zdcMaxAmp;
    bool cXn = coll.energyCommonZNC() > zdcMaxAmp && coll.timeZNC() < zdcMaxTime;
    bool c0n = coll.energyCommonZNC() <= zdcMaxAmp;

    if (this->neutronClass.value == "XnXn") {
      if (aXn && cXn) {
        return true;
      } else {
        return false;
      }
    } else if (this->neutronClass.value == "Xn0n") {
      if (aXn && c0n) {
        return true;
      } else {
        return false;
      }
    } else if (this->neutronClass.value == "0nXn") {
      if (a0n && cXn) {
        return true;
      } else {
        return false;
      }
    } else if (this->neutronClass.value == "0n0n") {
      if (a0n && c0n) {
        return true;
      } else {
        return false;
      }
    } else {
      // "Any" class
      return true;
    }
  } // End of Neutron class selection function

  template <typename T>
  bool isSelectedTrack(T const& track,
                       float ptcut,
                       float etaCut,
                       float dcaxycut,
                       float dcazcut,
                       bool ifITS,
                       bool ifTPC,
                       float itschi2nclscut,
                       float tpcchi2nclscut,
                       float tpcNClsCrossedRowscut)
  {
    PxPyPzMVector trackVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);
    // pt cut
    if (trackVector.Pt() < ptcut) {
      return false;
    }
    // eta cut
    if (std::fabs(trackVector.Eta()) > etaCut) {
      return false;
    }
    // DCA Z cut
    if (std::fabs(track.dcaZ()) > dcazcut) {
      return false;
    }
    // DCA XY cut
    float maxDCAxy = 0.0105 + 0.035 / std::pow(trackVector.Pt(), 1.1);
    if (dcaxycut == 0 && (std::fabs(track.dcaXY()) > maxDCAxy)) {
      return false;
    } else if (dcaxycut != 0 && (std::fabs(track.dcaXY()) > dcaxycut)) {
      return false;
    }
    // ITS Track only
    if (ifITS && !track.hasITS()) {
      return false;
    }
    // TPC Track only
    if (ifTPC && !track.hasTPC()) {
      return false;
    }
    // ITS Chi2 per N Clusters cut
    if (track.hasITS() && track.itsChi2NCl() > itschi2nclscut) {
      return false;
    }
    // TPC Chi2 N Clusters cut
    if (track.hasTPC() && track.tpcChi2NCl() > tpcchi2nclscut) {
      return false;
    }
    // TPC N Clusters Findable cut
    if (track.hasTPC() && track.tpcNClsCrossedRows() < tpcNClsCrossedRowscut) {
      return false;
    }
    if (useOnlyPVtracks && !track.isPVContributor()) {
      return false;
    }
    // All cuts passed
    return true;
  } // End of Track Selection function

  template <typename T>
  bool ifPion(const T& candidate, bool use_tof, float nsigmatpc_cut, float nsigmatof_cut, bool ifCircularNSigmaCut)
  {
    if (ifCircularNSigmaCut) {
      if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmatof_cut * nsigmatof_cut)) {
        return true;
      }

      if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
        return true;
      }

      if (!use_tof && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
        return true;
      }
      return false;
    } else {
      if (use_tof && candidate.hasTOF() && (std::abs(candidate.tofNSigmaPi()) < nsigmatpc_cut) && (std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut)) {
        return true;
      }

      if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
        return true;
      }

      if (!use_tof && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
        return true;
      }
      return false;
    }
  }

  int getRunNumberIndex(int runNumber)
  {
    for (int i = 0; i < numRunNums; ++i) {
      if (runNos[i] == runNumber) {
        return i;
      }
    }
    return -1; // Not found
  } // End of getRunNumberIndex function

  std::string strFormat(double value, int precision = 2)
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
  }

  void setHistBinLabels()
  {

    // Event cuts labels
    std::string eventLabels[13] = {
      "No Cuts",
      "isCBTHadronOk",
      "UPC or STD",
      "vtxITSTPC=" + strFormat(vtxITSTPCcut, 0),
      "sbp=" + strFormat(sbpCut, 0),
      "itsROFb=" + strFormat(itsROFbCut, 0),
      "tfb=" + strFormat(tfbCut, 0),
      "FT0A<=" + strFormat(fv0Cut),
      "FT0C<=" + strFormat(ft0cCut),
      "FV0A<=" + strFormat(ft0aCut),
      "Neutron Class: " + neutronClass.value,
      "n PV Contrib = 4",
      "V_{z} < " + strFormat(vZCut) + " cm"};
    int numEventCuts = 13;

    // Tracks cuts labels
    std::string trackLabels[12] = {
      "No Cuts",
      "pT>" + strFormat(pTcut) + " GeV/c",
      "|#eta|<" + strFormat(etaCut),
      "DCA Z<" + strFormat(dcaZcut) + " cm",
      "DCA XY cut",
      "hasITS = " + std::to_string(useITStracksOnly),
      "hasTPC = " + std::to_string(useTPCtracksOnly),
      "itsChi2NCl<" + strFormat(itsChi2NClsCut),
      "tpcChi2NCl<" + strFormat(tpcChi2NClsCut),
      "tpcNClsCrossedRows>" + strFormat(tpcNClsCrossedRowsCut),
      "#pi tracks (TPC+TOF)",
      "isPVContributor"};
    int numTrackCuts = 12;

    auto h1 = histosDataCounter.get<TH2>(HIST("EventsCounts_vs_runNo"));
    auto h2 = histosDataCounter.get<TH2>(HIST("TracksCounts_vs_runNo"));
    auto h3 = histos4piKin.get<THnSparse>(HIST("zero-charge"));
    auto h4 = histos4piKin.get<THnSparse>(HIST("non-zero-charge"));
    auto h5 = histosMCtruth.get<THnSparse>(HIST("Four-pion"));

    for (int i = 0; i < numEventCuts; ++i) {
      h1->GetYaxis()->SetBinLabel(i + 1, eventLabels[i].c_str());
    }
    for (int i = 0; i < numTrackCuts; ++i) {
      h2->GetYaxis()->SetBinLabel(i + 1, trackLabels[i].c_str());
    }
    for (int i = 0; i < numRunNums; ++i) {
      std::string runLabel = std::to_string(runNos[i]);
      h1->GetXaxis()->SetBinLabel(i + 1, runLabel.c_str());
      h2->GetXaxis()->SetBinLabel(i + 1, runLabel.c_str());
      h3->GetAxis(7)->SetBinLabel(i + 1, runLabel.c_str());
      h4->GetAxis(5)->SetBinLabel(i + 1, runLabel.c_str());
      h5->GetAxis(5)->SetBinLabel(i + 1, runLabel.c_str());
    }

  } // end of setHistBinLabels function

}; // End of Struct exclusiveRhoTo4Pi

int ExclusiveRhoTo4Pi::runNos[113] = {
  544013, 544028, 544032, 544091, 544095, 544098, 544116, 544121, 544122, 544123,
  544124, 544184, 544185, 544389, 544390, 544391, 544392, 544451, 544454, 544474,
  544475, 544476, 544477, 544490, 544491, 544492, 544508, 544510, 544511, 544512,
  544514, 544515, 544518, 544548, 544549, 544550, 544551, 544564, 544565, 544567,
  544568, 544580, 544582, 544583, 544585, 544614, 544640, 544652, 544653, 544672,
  544674, 544692, 544693, 544694, 544696, 544739, 544742, 544754, 544767, 544794,
  544795, 544797, 544813, 544868, 544886, 544887, 544896, 544911, 544913, 544914,
  544917, 544931, 544947, 544961, 544963, 544964, 544968, 544991, 544992, 545004,
  545008, 545009, 545041, 545042, 545044, 545047, 545060, 545062, 545063, 545064,
  545066, 545086, 545103, 545117, 545171, 545184, 545185, 545210, 545222, 545223,
  545246, 545249, 545262, 545289, 545291, 545294, 545295, 545296, 545311, 545312,
  545332, 545345, 545367};

int ExclusiveRhoTo4Pi::numRunNums = 113;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusiveRhoTo4Pi>(cfgc)};
}
