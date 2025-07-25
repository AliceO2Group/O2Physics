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

#include "Common/DataModel/PIDResponse.h"

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
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace branch
{
// Run Number
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
// Check UPC mode
DECLARE_SOA_COLUMN(IfCheckUPCmode, ifCheckUPCmode, uint16_t);
// vertex Position
DECLARE_SOA_COLUMN(PosX, posX, double);
DECLARE_SOA_COLUMN(PosY, posY, double);
DECLARE_SOA_COLUMN(PosZ, posZ, double);
// FIT signals
DECLARE_SOA_COLUMN(Fv0signal, fv0signal, double);
DECLARE_SOA_COLUMN(Ft0asignal, ft0asignal, double);
DECLARE_SOA_COLUMN(Ft0csignal, ft0csignal, double);
DECLARE_SOA_COLUMN(Fddasignal, fddasignal, double);
DECLARE_SOA_COLUMN(Fddcsignal, fddcsignal, double);
// FIT times
DECLARE_SOA_COLUMN(TimeFv0, timeFv0, double);
DECLARE_SOA_COLUMN(TimeFt0a, timeFt0a, double);
DECLARE_SOA_COLUMN(TimeFt0c, timeFt0c, double);
DECLARE_SOA_COLUMN(TimeFdda, timeFdda, double);
DECLARE_SOA_COLUMN(TimeFddc, timeFddc, double);
// ZDC times
DECLARE_SOA_COLUMN(TimeZna, timeZna, double);
DECLARE_SOA_COLUMN(TimeZnc, timeZnc, double);
// Occupancy
DECLARE_SOA_COLUMN(Occupancy, occupancy, double);
// Atleast one pion has TOF
DECLARE_SOA_COLUMN(HasAtLeastOneTOF, hasAtLeastOneTOF, bool);
// DCA XY
DECLARE_SOA_COLUMN(Dcaxy1, dcaxy1, double);
DECLARE_SOA_COLUMN(Dcaxy2, dcaxy2, double);
DECLARE_SOA_COLUMN(Dcaxy3, dcaxy3, double);
DECLARE_SOA_COLUMN(Dcaxy4, dcaxy4, double);
// DCA Z
DECLARE_SOA_COLUMN(Dcaz1, dcaz1, double);
DECLARE_SOA_COLUMN(Dcaz2, dcaz2, double);
DECLARE_SOA_COLUMN(Dcaz3, dcaz3, double);
DECLARE_SOA_COLUMN(Dcaz4, dcaz4, double);
// TPC nSigmaPi
DECLARE_SOA_COLUMN(TpcNsigmaPi1, tpcNsigmaPi1, double);
DECLARE_SOA_COLUMN(TpcNsigmaPi2, tpcNsigmaPi2, double);
DECLARE_SOA_COLUMN(TpcNsigmaPi3, tpcNsigmaPi3, double);
DECLARE_SOA_COLUMN(TpcNsigmaPi4, tpcNsigmaPi4, double);
// TPC nSigmaKa
DECLARE_SOA_COLUMN(TpcNsigmaKa1, tpcNsigmaKa1, double);
DECLARE_SOA_COLUMN(TpcNsigmaKa2, tpcNsigmaKa2, double);
DECLARE_SOA_COLUMN(TpcNsigmaKa3, tpcNsigmaKa3, double);
DECLARE_SOA_COLUMN(TpcNsigmaKa4, tpcNsigmaKa4, double);
// TPC nSigmaPr
DECLARE_SOA_COLUMN(TpcNsigmaPr1, tpcNsigmaPr1, double);
DECLARE_SOA_COLUMN(TpcNsigmaPr2, tpcNsigmaPr2, double);
DECLARE_SOA_COLUMN(TpcNsigmaPr3, tpcNsigmaPr3, double);
DECLARE_SOA_COLUMN(TpcNsigmaPr4, tpcNsigmaPr4, double);
// TPC nSigmaEl
DECLARE_SOA_COLUMN(TpcNsigmaEl1, tpcNsigmaEl1, double);
DECLARE_SOA_COLUMN(TpcNsigmaEl2, tpcNsigmaEl2, double);
DECLARE_SOA_COLUMN(TpcNsigmaEl3, tpcNsigmaEl3, double);
DECLARE_SOA_COLUMN(TpcNsigmaEl4, tpcNsigmaEl4, double);
// TPC nSigmaMu
DECLARE_SOA_COLUMN(TpcNsigmaMu1, tpcNsigmaMu1, double);
DECLARE_SOA_COLUMN(TpcNsigmaMu2, tpcNsigmaMu2, double);
DECLARE_SOA_COLUMN(TpcNsigmaMu3, tpcNsigmaMu3, double);
DECLARE_SOA_COLUMN(TpcNsigmaMu4, tpcNsigmaMu4, double);
// TPC Chi2
DECLARE_SOA_COLUMN(TpcChi21, tpcChi21, double);
DECLARE_SOA_COLUMN(TpcChi22, tpcChi22, double);
DECLARE_SOA_COLUMN(TpcChi23, tpcChi23, double);
DECLARE_SOA_COLUMN(TpcChi24, tpcChi24, double);
// TPC NClsFindable
DECLARE_SOA_COLUMN(TpcNClsFindable1, tpcNClsFindable1, double);
DECLARE_SOA_COLUMN(TpcNClsFindable2, tpcNClsFindable2, double);
DECLARE_SOA_COLUMN(TpcNClsFindable3, tpcNClsFindable3, double);
DECLARE_SOA_COLUMN(TpcNClsFindable4, tpcNClsFindable4, double);
// ITS Chi2
DECLARE_SOA_COLUMN(ItsChi21, itsChi21, double);
DECLARE_SOA_COLUMN(ItsChi22, itsChi22, double);
DECLARE_SOA_COLUMN(ItsChi23, itsChi23, double);
DECLARE_SOA_COLUMN(ItsChi24, itsChi24, double);
// PionPt
DECLARE_SOA_COLUMN(PionPt1, pionPt1, double);
DECLARE_SOA_COLUMN(PionPt2, pionPt2, double);
DECLARE_SOA_COLUMN(PionPt3, pionPt3, double);
DECLARE_SOA_COLUMN(PionPt4, pionPt4, double);
// Pion Eta
DECLARE_SOA_COLUMN(PionEta1, pionEta1, double);
DECLARE_SOA_COLUMN(PionEta2, pionEta2, double);
DECLARE_SOA_COLUMN(PionEta3, pionEta3, double);
DECLARE_SOA_COLUMN(PionEta4, pionEta4, double);
// Pion Phi
DECLARE_SOA_COLUMN(PionPhi1, pionPhi1, double);
DECLARE_SOA_COLUMN(PionPhi2, pionPhi2, double);
DECLARE_SOA_COLUMN(PionPhi3, pionPhi3, double);
DECLARE_SOA_COLUMN(PionPhi4, pionPhi4, double);
// Pion Rapidity
DECLARE_SOA_COLUMN(PionRapidity1, pionRapidity1, double);
DECLARE_SOA_COLUMN(PionRapidity2, pionRapidity2, double);
DECLARE_SOA_COLUMN(PionRapidity3, pionRapidity3, double);
DECLARE_SOA_COLUMN(PionRapidity4, pionRapidity4, double);
// Four Pion Pt, Eta, Phi Rapidity
DECLARE_SOA_COLUMN(FourPionPt, fourPionPt, double);
DECLARE_SOA_COLUMN(FourPionEta, fourPionEta, double);
DECLARE_SOA_COLUMN(FourPionPhi, fourPionPhi, double);
DECLARE_SOA_COLUMN(FourPionRapidity, fourPionRapidity, double);
DECLARE_SOA_COLUMN(FourPionMass, fourPionMass, double);
// Four Pion Phi Pair 1, Pair 2, CosTheta Pair 1, CosTheta Pair 2
DECLARE_SOA_COLUMN(FourPionPhiPair1, fourPionPhiPair1, double);
DECLARE_SOA_COLUMN(FourPionPhiPair2, fourPionPhiPair2, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair1, fourPionCosThetaPair1, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair2, fourPionCosThetaPair2, double);
} // namespace branch

DECLARE_SOA_TABLE(SignalData, "AOD", "signalData",
                  branch::RunNumber,

                  branch::IfCheckUPCmode,

                  branch::PosX,
                  branch::PosY,
                  branch::PosZ,

                  branch::Fv0signal,
                  branch::Ft0asignal,
                  branch::Ft0csignal,
                  branch::Fddasignal,
                  branch::Fddcsignal,

                  branch::TimeFv0,
                  branch::TimeFt0a,
                  branch::TimeFt0c,
                  branch::TimeFdda,
                  branch::TimeFddc,
                  branch::TimeZna,
                  branch::TimeZnc,

                  branch::Occupancy,

                  branch::HasAtLeastOneTOF,

                  branch::Dcaxy1,
                  branch::Dcaxy2,
                  branch::Dcaxy3,
                  branch::Dcaxy4,

                  branch::Dcaz1,
                  branch::Dcaz2,
                  branch::Dcaz3,
                  branch::Dcaz4,

                  branch::TpcNsigmaPi1,
                  branch::TpcNsigmaPi2,
                  branch::TpcNsigmaPi3,
                  branch::TpcNsigmaPi4,

                  branch::TpcNsigmaKa1,
                  branch::TpcNsigmaKa2,
                  branch::TpcNsigmaKa3,
                  branch::TpcNsigmaKa4,

                  branch::TpcNsigmaPr1,
                  branch::TpcNsigmaPr2,
                  branch::TpcNsigmaPr3,
                  branch::TpcNsigmaPr4,

                  branch::TpcNsigmaEl1,
                  branch::TpcNsigmaEl2,
                  branch::TpcNsigmaEl3,
                  branch::TpcNsigmaEl4,

                  branch::TpcNsigmaMu1,
                  branch::TpcNsigmaMu2,
                  branch::TpcNsigmaMu3,
                  branch::TpcNsigmaMu4,

                  branch::TpcChi21,
                  branch::TpcChi22,
                  branch::TpcChi23,
                  branch::TpcChi24,

                  branch::TpcNClsFindable1,
                  branch::TpcNClsFindable2,
                  branch::TpcNClsFindable3,
                  branch::TpcNClsFindable4,

                  branch::ItsChi21,
                  branch::ItsChi22,
                  branch::ItsChi23,
                  branch::ItsChi24,

                  branch::PionPt1,
                  branch::PionPt2,
                  branch::PionPt3,
                  branch::PionPt4,

                  branch::PionEta1,
                  branch::PionEta2,
                  branch::PionEta3,
                  branch::PionEta4,

                  branch::PionPhi1,
                  branch::PionPhi2,
                  branch::PionPhi3,
                  branch::PionPhi4,

                  branch::PionRapidity1,
                  branch::PionRapidity2,
                  branch::PionRapidity3,
                  branch::PionRapidity4,

                  branch::FourPionPt,
                  branch::FourPionEta,
                  branch::FourPionPhi,
                  branch::FourPionRapidity,
                  branch::FourPionMass,
                  branch::FourPionPhiPair1,
                  branch::FourPionPhiPair2,
                  branch::FourPionCosThetaPair1,
                  branch::FourPionCosThetaPair2);

DECLARE_SOA_TABLE(BkgroundData, "AOD", "bkgroundData",
                  branch::RunNumber,

                  branch::IfCheckUPCmode,

                  branch::PosX,
                  branch::PosY,
                  branch::PosZ,

                  branch::Fv0signal,
                  branch::Ft0asignal,
                  branch::Ft0csignal,
                  branch::Fddasignal,
                  branch::Fddcsignal,

                  branch::TimeFv0,
                  branch::TimeFt0a,
                  branch::TimeFt0c,
                  branch::TimeFdda,
                  branch::TimeFddc,
                  branch::TimeZna,
                  branch::TimeZnc,

                  branch::Occupancy,

                  branch::HasAtLeastOneTOF,

                  branch::Dcaxy1,
                  branch::Dcaxy2,
                  branch::Dcaxy3,
                  branch::Dcaxy4,

                  branch::Dcaz1,
                  branch::Dcaz2,
                  branch::Dcaz3,
                  branch::Dcaz4,

                  branch::TpcNsigmaPi1,
                  branch::TpcNsigmaPi2,
                  branch::TpcNsigmaPi3,
                  branch::TpcNsigmaPi4,

                  branch::TpcNsigmaKa1,
                  branch::TpcNsigmaKa2,
                  branch::TpcNsigmaKa3,
                  branch::TpcNsigmaKa4,

                  branch::TpcNsigmaPr1,
                  branch::TpcNsigmaPr2,
                  branch::TpcNsigmaPr3,
                  branch::TpcNsigmaPr4,

                  branch::TpcNsigmaEl1,
                  branch::TpcNsigmaEl2,
                  branch::TpcNsigmaEl3,
                  branch::TpcNsigmaEl4,

                  branch::TpcNsigmaMu1,
                  branch::TpcNsigmaMu2,
                  branch::TpcNsigmaMu3,
                  branch::TpcNsigmaMu4,

                  branch::TpcChi21,
                  branch::TpcChi22,
                  branch::TpcChi23,
                  branch::TpcChi24,

                  branch::TpcNClsFindable1,
                  branch::TpcNClsFindable2,
                  branch::TpcNClsFindable3,
                  branch::TpcNClsFindable4,

                  branch::ItsChi21,
                  branch::ItsChi22,
                  branch::ItsChi23,
                  branch::ItsChi24,

                  branch::PionPt1,
                  branch::PionPt2,
                  branch::PionPt3,
                  branch::PionPt4,

                  branch::PionEta1,
                  branch::PionEta2,
                  branch::PionEta3,
                  branch::PionEta4,

                  branch::PionPhi1,
                  branch::PionPhi2,
                  branch::PionPhi3,
                  branch::PionPhi4,

                  branch::PionRapidity1,
                  branch::PionRapidity2,
                  branch::PionRapidity3,
                  branch::PionRapidity4,

                  branch::FourPionPt,
                  branch::FourPionEta,
                  branch::FourPionPhi,
                  branch::FourPionRapidity,
                  branch::FourPionMass);
} // namespace o2::aod

struct ExclusiveRhoTo4Pi {
  SGSelector sgSelector;
  // Defining constants
  int numFourPionTracks = 4;
  int numPiPlus = 2;
  int numPiMinus = 2;
  float zeroPointEight = 0.8;
  // Derived Data
  Produces<aod::SignalData> sigFromData;
  Produces<aod::BkgroundData> bkgFromData;
  // Histogram Registry
  HistogramRegistry histosData{"histosData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable Event parameters
  Configurable<bool> ifCheckUPCmode{"ifCheckUPCmode", false, "Enable UPC reconstruction only"};
  Configurable<float> vZCut{"vZCut", 10., "Vertex Cut"};
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 150., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> zdcCut{"zdcCut", 0., "ZDC threshold"};
  Configurable<uint16_t> numPVContrib{"numPVContrib", 4, "Number of PV Contributors"};
  Configurable<int> gapSideCut{"gapSideCut", 2, "Gap Side"};
  Configurable<int> sbpCut{"sbpCut", 1, "Sbp"};
  Configurable<int> itsROFbCut{"itsROFbCut", 1, "itsROFbCut"};
  Configurable<int> vtxITSTPCcut{"vtxITSTPCcut", 1, "vtxITSTPCcut"};
  Configurable<int> tfbCut{"tfbCut", 1, "tfbCut"};
  // track Selection mode
  Configurable<int> trackSelectionMode{"trackSelectionMode", 0, "Different modes of track selection"};
  // Configurable Track parameters common to mode 0 and 1
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use Only PV tracks"};
  Configurable<bool> useITS{"useITS", true, "only use tracks with hit in ITS"};
  Configurable<bool> useTPC{"useTPC", true, "has TPC hit"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min TPC Findable Clusters"};
  Configurable<float> pTcut{"pTcut", 0.1, "Track Pt"};
  Configurable<float> dcaZcut{"dcaZcut", 1, "dcaZ cut"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  // Configurable Track parameters for mode 0 only
  Configurable<uint8_t> itsNClsCut{"itsNClsCut", 4, "Min No of itsNCls"};
  Configurable<uint8_t> itsClusterMapCut{"itsClusterMapCut", 1, "min no of ITS clusters in cluster map"};
  Configurable<float> itsChi2NClCut{"itsChi2NClCut", 3.0, "Max ITS Chi2/NCl"};
  Configurable<float> minFoundTPCclusters{"minFoundTPCclusters", 120, "Min TPC Findable Clusters"};
  Configurable<float> tpcChi2NClsMin{"tpcChi2NClsMin", 1.0, "Min TPC Chi2/NCls"};
  Configurable<float> tpcChi2NClsMax{"tpcChi2NClsMax", 3.0, "Max TPC Chi2/NCls"};
  Configurable<float> tpcNClsCrossedRowsCut{"tpcNClsCrossedRowsCut", 130, "Min TPC Crossed Rows"};
  Configurable<float> tpcCrossedRowsOverFindableCut{"tpcCrossedRowsOverFindableCut", 1.0, "Min TPC Crossed Rows over Findable Clusters"};
  // Configurable Track parameters for mode: 1 only
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "ITS Chi2"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4.0, "TPC Chi2"};
  // Configurable PID parameters
  Configurable<bool> useTOF{"useTOF", true, "has TOF for PID"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 3, "TPC cut"};
  Configurable<float> nSigmaTOFcut{"nSigmaTOFcut", 3, "TOF cut"};
  // Configurable Rho parameters
  Configurable<float> rhoRapCut{"rhoRapCut", 0.5, "Max abs Rapidity of rho"};
  Configurable<float> rhoPtCut{"rhoPtCut", 0.15, "Min Pt of rho"};
  Configurable<float> rhoMassMin{"rhoMassMin", 1, "Min Mass of rho"};
  Configurable<float> rhoMassMax{"rhoMassMax", 2.5, "Max Mass of rho"};
  // Axis Configurations
  ConfigurableAxis pTAxis{"pTAxis", {1000, 0, 2}, "Axis for pT histograms"};
  ConfigurableAxis etaAxis{"etaAxis", {1000, -1.1, 1.1}, "Axis for Eta histograms"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {1000, -2.5, 2.5}, "Axis for Rapidity histograms"};
  ConfigurableAxis invMassAxis{"invMassAxis", {1000, 1, 2.5}, "Axis for Phi histograms"};
  ConfigurableAxis phiAxis{"phiAxis", {360, -1 * o2::constants::math::PI, o2::constants::math::PI}, "Axis for Phi histograms"};
  ConfigurableAxis cosThetaAxis{"cosThetaAxis", {360, -1, 1}, "Axis for cos Theta histograms"};

  void init(InitContext const&)
  {
    // QA plots: Event Counter
    histosData.add("EventsCounts_vs_runNo", "Number of Selected 4-Pion Events per Run; Run Number; Number of Events", kTH2F, {{1355, 544013, 545367}, {21, -1, 20}});
    histosData.add("TracksCounts_vs_runNo", "Number of Selected Tracks per Run; Run Number; Number of Tracks", kTH2F, {{1355, 544013, 545367}, {20, 0, 20}});
    // QA plots: event selection
    histosData.add("FT0A", "T0A amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosData.add("FT0C", "T0C amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosData.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosData.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosData.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});
    histosData.add("vertexX", "Vertex X; Vertex X [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosData.add("vertexY", "Vertex Y; Vertex Y [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosData.add("vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{2000, -15, 15}});
    histosData.add("occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{20000, 0, 20000}});
    histosData.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histosData.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    // QA plots: tracks
    histosData.add("dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{5000, -1, 1}});
    histosData.add("dcaXY_pions", "dcaXY_pions; dcaXY of Pions [cm]; Counts", kTH1F, {{5000, -1, 1}});
    histosData.add("dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{5000, -1, 1}});
    histosData.add("dcaZ_pions", "dcaZ_pions; dcaZ of Pions [cm]; Counts", kTH1F, {{5000, -1, 1}});
    histosData.add("tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("tpcNClsFindable", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    histosData.add("itsClusterMap", "ITS Cluster Map; itsClusterMap; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: PID
    histosData.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosData.add("tpcSignal_pions", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosData.add("tpcNSigmaPi_all", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPi_pions", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaKa_pions", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPr_pions", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaEl_pions", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaMu_pions", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histosData.add("tofBeta_pions", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histosData.add("tofNSigmaPi_all", "TOF nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaPi_pions", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaKa_pions", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaPr_pions", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaEl_pions", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaMu_pions", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    // Track Transverse Momentum
    histosData.add("pT_track_all", "pT with track selection; pT [GeV/c]; Counts", kTH1F, {pTAxis});
    histosData.add("pT_track_pions", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosData.add("pT_track_pions_contributed", "pT with track selection and PID selection of Pi which are contributed to selected event; pT [GeV/c]; Events", kTH1F, {pTAxis});
    // Track Pseudorapidity
    histosData.add("eta_track_all", "Pseudorapidity with track selection; #eta; Counts", kTH1F, {etaAxis});
    histosData.add("eta_track_pions", "Pseudorapidity with track selection and PID selection of Pi; #eta; Events", kTH1F, {etaAxis});
    histosData.add("eta_track_pions_contributed", "Pseudorapidity with track selection and PID selection of Pi which are contributed to selected event; #eta; Events", kTH1F, {etaAxis});
    // Track Phi
    histosData.add("phi_track_all", "Phi with track selection; #phi [rad]; Counts", kTH1F, {phiAxis});
    histosData.add("phi_track_pions", "Phi with track selection and PID selection of Pi; #phi [rad]; Events", kTH1F, {phiAxis});
    histosData.add("phi_track_pions_contributed", "Phi with track selection and PID selection of Pi which are contributed to selected event; #phi [rad]; Events", kTH1F, {phiAxis});
    // Track Rapidity
    histosData.add("rapidity_track_all", "Rapidity with track selection; y; Counts", kTH1F, {rapidityAxis});
    histosData.add("rapidity_track_pions", "Rapidity with track selection and PID selection of Pi; y; Events", kTH1F, {rapidityAxis});
    histosData.add("rapidity_track_pions_contributed", "Rapidity with track selection and PID selection of Pi which are contributed to selected event; y; Events", kTH1F, {rapidityAxis});
    // Four Pion Transverse Momentum
    histosData.add("fourpion_pT_0_charge", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosData.add("fourpion_pT_0_charge_within_rap", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosData.add("fourpion_pT_non_0_charge", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosData.add("fourpion_pT_non_0_charge_within_rap", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    // Four Pion Eta
    histosData.add("fourpion_eta_0_charge", "Four Pion #eta (0 charge); #eta; Events", kTH1F, {etaAxis});
    histosData.add("fourpion_eta_0_charge_within_rap", "Four Pion #eta (0 charge within rap); #eta; Events", kTH1F, {etaAxis});
    histosData.add("fourpion_eta_non_0_charge", "Four Pion #eta (non 0 charge); #eta; #eta; Events", kTH1F, {etaAxis});
    histosData.add("fourpion_eta_non_0_charge_within_rap", "Four Pion #eta (non 0 charge within rap); #eta; Events", kTH1F, {etaAxis});
    // Four Pion Phi
    histosData.add("fourpion_phi_0_charge", "Four Pion #phi (0 charge); #phi [rad]; Events", kTH1F, {phiAxis});
    histosData.add("fourpion_phi_0_charge_within_rap", "Four Pion #phi (0 charge within rap); #phi [rad]; Events", kTH1F, {phiAxis});
    histosData.add("fourpion_phi_non_0_charge", "Four Pion #phi (non 0 charge); #phi [rad]; Events", kTH1F, {phiAxis});
    histosData.add("fourpion_phi_non_0_charge_within_rap", "Four Pion #phi (non 0 charge within rap); #phi [rad]; Events", kTH1F, {phiAxis});
    // Four Pion Rapidity
    histosData.add("fourpion_rap_0_charge", "Four Pion Rapidity (0 charge); y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosData.add("fourpion_rap_0_charge_within_rap", "Four Pion Rapidity (0 charge within rap); y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosData.add("fourpion_rap_non_0_charge", "Four Pion Rapidity (non 0 charge); y; Events", kTH1F, {rapidityAxis});
    histosData.add("fourpion_rap_non_0_charge_within_rap", "Four Pion Rapidity (non 0 charge within rap); y; Events", kTH1F, {rapidityAxis});
    // Four Pion Mass
    histosData.add("fourpion_mass_0_charge", "Four Pion Invariant Mass (0 charge); m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]; Events", kTH1F, {invMassAxis});
    histosData.add("fourpion_mass_0_charge_within_rap", "Four Pion Invariant Mass (0 charge within rap); m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]; Events", kTH1F, {invMassAxis});
    histosData.add("fourpion_mass_non_0_charge", "Four Pion Invariant Mass (non 0 charge); m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]; Events", kTH1F, {invMassAxis});
    histosData.add("fourpion_mass_non_0_charge_within_rap", "Four Pion Invariant Mass (non 0 charge within rap); m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]; Events", kTH1F, {invMassAxis});
    // Pair Invariant Mass
    histosData.add("twopion_mass_1", "Invariant Mass Distribution of 2 pions 1 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_2", "Invariant Mass Distribution of 2 pions 2 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_3", "Invariant Mass Distribution of 2 pions 3 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_4", "Invariant Mass Distribution of 2 pions 4 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    // Four Pion Invariant Mass
    histosData.add("fourpion_mass_0_charge_domA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});               // pT < 0.15GeV
    histosData.add("fourpion_mass_0_charge_domB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});         // 0.15GeV < pT < 0.8GeV
    histosData.add("fourpion_mass_0_charge_domC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});               // 0.8GeV < pT
    histosData.add("fourpion_mass_non_0_charge_domA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // pT < 0.15GeV
    histosData.add("fourpion_mass_non_0_charge_domB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis}); // 0.15GeV < pT < 0.8GeV
    histosData.add("fourpion_mass_non_0_charge_domC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // 0.8GeV < pT
    // Collin Soper Theta and Phi
    histosData.add("collin_soper_phi_1", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_phi_2", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_costheta_1", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("collin_soper_costheta_2", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("phi_vs_costheta_1", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    histosData.add("phi_vs_costheta_2", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    // Collin Soper Theta and Phi after selection
    histosData.add("collin_soper_phi_small_mass", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_phi_large_mass", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_costheta_small_mass", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("collin_soper_costheta_large_mass", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("phi_vs_costheta_small_mass", "Phi vs cosTheta for small mass; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    histosData.add("phi_vs_costheta_large_mass", "Phi vs cosTheta for large mass; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
  } // End of init function

  using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;

  void processData(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {

    // no cuts
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), -1);

    // Check if the event is in UPC mode
    if (ifCheckUPCmode && (collision.flags() != 1)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 0);

    // FTOA
    if (!(collision.totalFT0AmplitudeA() <= ft0aCut)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 1);

    // FT0C
    if (!(collision.totalFT0AmplitudeC() <= ft0cCut)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 2);

    // FV0
    if (!(collision.totalFV0AmplitudeA() <= fv0Cut)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 3);

    // noSamebunchPileup
    if (collision.sbp() != sbpCut) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 4);

    // kIsVertexITSTPC
    if (collision.vtxITSTPC() != vtxITSTPCcut) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 5);

    // kNoITSROFrameBorder
    if (collision.itsROFb() != itsROFbCut) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 6);

    // kNoTimeFrameBorder
    if (collision.tfb() != tfbCut) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 7);

    // ZDC
    if (!(collision.energyCommonZNA() <= zdcCut || collision.energyCommonZNC() <= zdcCut)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 8);

    // Vertex Z cut
    if (!(std::abs(collision.posZ()) <= vZCut)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 9);

    // true Gap Side
    if (sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut) != gapSideCut) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 10);

    // number of PV contributors
    if (!(collision.numContrib() == numPVContrib)) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 11);

    histosData.fill(HIST("GapSide"), collision.gapSide());
    histosData.fill(HIST("TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
    histosData.fill(HIST("vertexX"), collision.posX());
    histosData.fill(HIST("vertexY"), collision.posY());
    histosData.fill(HIST("vertexZ"), collision.posZ());
    histosData.fill(HIST("occupancy"), collision.occupancyInTime());
    histosData.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histosData.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histosData.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histosData.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histosData.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    std::vector<decltype(tracks.begin())> selectedTracks;
    std::vector<decltype(tracks.begin())> selectedPionTracks;
    std::vector<decltype(tracks.begin())> selectedPionPlusTracks;
    std::vector<decltype(tracks.begin())> selectedPionMinusTracks;

    for (const auto& t0 : tracks) {

      ROOT::Math::PxPyPzMVector trackVector(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
      // no Cuts
      histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 0);

      if (trackSelectionMode == 0) {

        // is PV Contributor
        if (!(t0.isPVContributor() == useOnlyPVtracks)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 1);

        // has ITS hit
        if ((useITS == true) && (t0.hasITS() != true)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 2);

        // min no of itsNCls
        if (t0.itsNCls() < itsNClsCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 3);

        // min ITS chi2NCl
        if (t0.itsChi2NCl() > itsChi2NClCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 4);

        // has TPC hit
        if ((useTPC == true) && (t0.hasTPC() != true)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 5);

        // min no of found TPC clusters
        if (t0.tpcNClsFindable() - t0.tpcNClsFindableMinusFound() < minFoundTPCclusters) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 6);

        // range of tpcChi2NCl
        if (!((tpcChi2NClsMin < t0.tpcChi2NCl()) && (t0.tpcChi2NCl() < tpcChi2NClsMax))) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 7);

        // tpcNClsCrossedRows
        if (t0.tpcNClsCrossedRows() < tpcNClsCrossedRowsCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 8);

        // ratio of crossed TPC rows over findable clusters
        if ((t0.tpcNClsCrossedRows() / t0.tpcNClsFindable()) < tpcCrossedRowsOverFindableCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 9);

        // pT cut
        if (trackVector.Pt() < pTcut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 10);

        // dcaZ cut
        if ((std::abs(t0.dcaZ()) > dcaZcut) || (t0.dcaXY() > getMaxDCAxy(trackVector.Pt()))) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 11);

        // eta cut
        if (std::abs(trackVector.Eta()) > etaCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 12);
      } // end of trackSelectionMode == 0

      if (trackSelectionMode == 1) {
        // is PV Contributor
        if (!(t0.isPVContributor() == useOnlyPVtracks)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 1);

        // pT cut
        if (trackVector.Pt() < pTcut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 2);

        // eta cut
        if (std::abs(trackVector.Eta()) > etaCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 3);

        // dcaZ cut
        if ((std::abs(t0.dcaZ()) > dcaZcut)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 4);

        // dcaXY cut
        if (std::abs(t0.dcaXY()) > getMaxDCAxy(trackVector.Pt())) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 5);

        // has ITS hit
        if ((useITS == true) && (t0.hasITS() != true)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 6);

        // has TPC hit
        if ((useTPC == true) && (t0.hasTPC() != true)) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 7);

        // ITS Chi2 Cut
        if (t0.itsChi2NCl() > itsChi2Cut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 8);

        // TPC Chi2 Cut
        if (t0.tpcChi2NCl() > tpcChi2Cut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 9);

        // TPC Clusters findable cut
        if (t0.tpcNClsFindable() < tpcNClsFindableCut) {
          continue;
        }
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 10);
      } // end of trackSelectionMode == 1

      selectedTracks.push_back(t0);
      if (selectionPIDPion(t0, useTOF, nSigmaTPCcut, nSigmaTOFcut)) {
        selectedPionTracks.push_back(t0);
        histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 13);
        if (t0.sign() == 1) {
          selectedPionPlusTracks.push_back(t0);
          histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 14);
        }
        if (t0.sign() == -1) {
          selectedPionMinusTracks.push_back(t0);
          histosData.fill(HIST("TracksCounts_vs_runNo"), collision.runNumber(), 15);
        }
      } // End of Selection PID Pion
    } // End of loop over tracks

    int numSelectedTracks = static_cast<int>(selectedTracks.size());
    int numSelectedPionTracks = static_cast<int>(selectedPionTracks.size());
    int numPiPlusTracks = static_cast<int>(selectedPionPlusTracks.size());
    int numPionMinusTracks = static_cast<int>(selectedPionMinusTracks.size());

    for (int i = 0; i < numSelectedTracks; i++) {
      ROOT::Math::PxPyPzMVector selectedTrackVector(selectedTracks[i].px(), selectedTracks[i].py(), selectedTracks[i].pz(), o2::constants::physics::MassPionCharged);
      histosData.fill(HIST("tpcSignal"), selectedTrackVector.P(), selectedTracks[i].tpcSignal());
      histosData.fill(HIST("tofBeta"), selectedTrackVector.P(), selectedTracks[i].beta());
      histosData.fill(HIST("tpcNSigmaPi_all"), selectedTracks[i].tpcNSigmaPi(), selectedTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaPi_all"), selectedTracks[i].tofNSigmaPi(), selectedTrackVector.Pt());
      histosData.fill(HIST("pT_track_all"), selectedTrackVector.Pt());
      histosData.fill(HIST("eta_track_all"), selectedTrackVector.Eta());
      histosData.fill(HIST("phi_track_all"), selectedTrackVector.Phi());
      histosData.fill(HIST("rapidity_track_all"), selectedTrackVector.Rapidity());
      histosData.fill(HIST("itsChi2NCl"), selectedTracks[i].itsChi2NCl());
      histosData.fill(HIST("tpcChi2NCl"), selectedTracks[i].tpcChi2NCl());
      histosData.fill(HIST("tpcNClsFindable"), selectedTracks[i].tpcNClsFindable());
      histosData.fill(HIST("dcaXY"), selectedTracks[i].dcaXY());
      histosData.fill(HIST("dcaZ"), selectedTracks[i].dcaZ());
      histosData.fill(HIST("itsClusterMap"), selectedTracks[i].itsClusterMap());
    } // End of loop over tracks with selection only

    for (int i = 0; i < numSelectedPionTracks; i++) {
      ROOT::Math::PxPyPzMVector selectedPionTrackVector(selectedPionTracks[i].px(), selectedPionTracks[i].py(), selectedPionTracks[i].pz(), o2::constants::physics::MassPionCharged);
      histosData.fill(HIST("tpcSignal_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].tpcSignal());
      histosData.fill(HIST("tofBeta_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].beta());
      histosData.fill(HIST("tpcNSigmaPi_pions"), selectedPionTracks[i].tpcNSigmaPi(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaKa_pions"), selectedPionTracks[i].tpcNSigmaKa(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaPr_pions"), selectedPionTracks[i].tpcNSigmaPr(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaEl_pions"), selectedPionTracks[i].tpcNSigmaEl(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaMu_pions"), selectedPionTracks[i].tpcNSigmaMu(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaPi_pions"), selectedPionTracks[i].tofNSigmaPi(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaKa_pions"), selectedPionTracks[i].tofNSigmaKa(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaPr_pions"), selectedPionTracks[i].tofNSigmaPr(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaEl_pions"), selectedPionTracks[i].tofNSigmaEl(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaMu_pions"), selectedPionTracks[i].tofNSigmaMu(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("pT_track_pions"), selectedPionTrackVector.Pt());
      histosData.fill(HIST("eta_track_pions"), selectedPionTrackVector.Eta());
      histosData.fill(HIST("phi_track_pions"), selectedPionTrackVector.Phi());
      histosData.fill(HIST("rapidity_track_pions"), selectedPionTrackVector.Rapidity());
      histosData.fill(HIST("dcaXY_pions"), selectedPionTracks[i].dcaXY());
      histosData.fill(HIST("dcaZ_pions"), selectedPionTracks[i].dcaZ());
    } // End of loop over tracks with selection and PID of pions

    // event should have exactly 4 pions
    if (numSelectedPionTracks != numFourPionTracks) {
      return;
    }
    histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 12);

    // Check if there is at least one track with TOF in the selected events (for derived Data)
    bool hasAtleastOneTOF = false;
    for (int i = 0; i < numPiPlusTracks; i++) {
      if (selectedPionPlusTracks[i].hasTOF() == true) {
        hasAtleastOneTOF = true;
        break;
      }
    }

    // Selecting Events with net charge = 0
    if (numPionMinusTracks == numPiMinus && numPiPlusTracks == numPiPlus) {

      histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 13);

      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      ROOT::Math::PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

      histosData.fill(HIST("pT_track_pions_contributed"), p1.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p2.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p3.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p4.Pt());

      histosData.fill(HIST("eta_track_pions_contributed"), p1.Eta());
      histosData.fill(HIST("eta_track_pions_contributed"), p2.Eta());
      histosData.fill(HIST("eta_track_pions_contributed"), p3.Eta());
      histosData.fill(HIST("eta_track_pions_contributed"), p4.Eta());

      histosData.fill(HIST("phi_track_pions_contributed"), p1.Phi());
      histosData.fill(HIST("phi_track_pions_contributed"), p2.Phi());
      histosData.fill(HIST("phi_track_pions_contributed"), p3.Phi());
      histosData.fill(HIST("phi_track_pions_contributed"), p4.Phi());

      histosData.fill(HIST("rapidity_track_pions_contributed"), p1.Rapidity());
      histosData.fill(HIST("rapidity_track_pions_contributed"), p2.Rapidity());
      histosData.fill(HIST("rapidity_track_pions_contributed"), p3.Rapidity());
      histosData.fill(HIST("rapidity_track_pions_contributed"), p4.Rapidity());

      k1.SetCoordinates(p1.Pt(), p1.Eta(), p1.Phi(), o2::constants::physics::MassPionCharged);
      k2.SetCoordinates(p2.Pt(), p2.Eta(), p2.Phi(), o2::constants::physics::MassPionCharged);
      k3.SetCoordinates(p3.Pt(), p3.Eta(), p3.Phi(), o2::constants::physics::MassPionCharged);
      k4.SetCoordinates(p4.Pt(), p4.Eta(), p4.Phi(), o2::constants::physics::MassPionCharged);

      ROOT::Math::PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
      k1234 = k1 + k2 + k3 + k4;

      k13 = k1 + k3;
      k14 = k1 + k4;
      k23 = k2 + k3;
      k24 = k2 + k4;

      histosData.fill(HIST("fourpion_pT_0_charge"), p1234.Pt());
      histosData.fill(HIST("fourpion_eta_0_charge"), p1234.Eta());
      histosData.fill(HIST("fourpion_phi_0_charge"), p1234.Phi());
      histosData.fill(HIST("fourpion_rap_0_charge"), p1234.Rapidity());
      histosData.fill(HIST("fourpion_mass_0_charge"), p1234.M());

      double fourPiPhiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
      double fourPiPhiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
      double fourPiCosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
      double fourPiCosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

      sigFromData(
        // run number
        collision.runNumber(),
        // UPC mode
        collision.flags(),
        // vertex
        collision.posX(), collision.posY(), collision.posZ(),
        // FIT Signals
        collision.totalFV0AmplitudeA(), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
        // FIT and ZDC Signals
        collision.timeFV0A(), collision.timeFT0A(), collision.timeFT0C(), collision.timeFDDA(), collision.timeFDDC(), collision.timeZNA(), collision.timeZNC(),
        // Occupancy
        collision.occupancyInTime(),
        // has atleast one TOF
        hasAtleastOneTOF,
        // DCA XY and Z
        selectedPionPlusTracks[0].dcaXY(), selectedPionPlusTracks[1].dcaXY(), selectedPionMinusTracks[0].dcaXY(), selectedPionMinusTracks[1].dcaXY(),
        selectedPionPlusTracks[0].dcaZ(), selectedPionPlusTracks[1].dcaZ(), selectedPionMinusTracks[0].dcaZ(), selectedPionMinusTracks[1].dcaZ(),
        // TPC N Sigma Pi
        selectedPionPlusTracks[0].tpcNSigmaPi(), selectedPionPlusTracks[1].tpcNSigmaPi(), selectedPionMinusTracks[0].tpcNSigmaPi(), selectedPionMinusTracks[1].tpcNSigmaPi(),
        // TPC N Sigma Ka
        selectedPionPlusTracks[0].tpcNSigmaKa(), selectedPionPlusTracks[1].tpcNSigmaKa(), selectedPionMinusTracks[0].tpcNSigmaKa(), selectedPionMinusTracks[1].tpcNSigmaKa(),
        // TPC N Sigma Pr
        selectedPionPlusTracks[0].tpcNSigmaPr(), selectedPionPlusTracks[1].tpcNSigmaPr(), selectedPionMinusTracks[0].tpcNSigmaPr(), selectedPionMinusTracks[1].tpcNSigmaPr(),
        // TPC N Sigma El
        selectedPionPlusTracks[0].tpcNSigmaEl(), selectedPionPlusTracks[1].tpcNSigmaEl(), selectedPionMinusTracks[0].tpcNSigmaEl(), selectedPionMinusTracks[1].tpcNSigmaEl(),
        // TPC N Sigma Mu
        selectedPionPlusTracks[0].tpcNSigmaMu(), selectedPionPlusTracks[1].tpcNSigmaMu(), selectedPionMinusTracks[0].tpcNSigmaMu(), selectedPionMinusTracks[1].tpcNSigmaMu(),
        // tpc Chi2 NCl
        selectedPionPlusTracks[0].tpcChi2NCl(), selectedPionPlusTracks[1].tpcChi2NCl(), selectedPionMinusTracks[0].tpcChi2NCl(), selectedPionMinusTracks[1].tpcChi2NCl(),
        // TPC NCls Findable
        selectedPionPlusTracks[0].tpcNClsFindable(), selectedPionPlusTracks[1].tpcNClsFindable(), selectedPionMinusTracks[0].tpcNClsFindable(), selectedPionMinusTracks[1].tpcNClsFindable(),
        // ITS Chi2 NCl
        selectedPionPlusTracks[0].itsChi2NCl(), selectedPionPlusTracks[1].itsChi2NCl(), selectedPionMinusTracks[0].itsChi2NCl(), selectedPionMinusTracks[1].itsChi2NCl(),
        // Pion Pt
        p1.Pt(), p2.Pt(), p3.Pt(), p4.Pt(),
        // Pion Eta
        p1.Eta(), p2.Eta(), p3.Eta(), p4.Eta(),
        // Pion Phi
        p1.Phi(), p2.Phi(), p3.Phi(), p4.Phi(),
        // Pion Rapidity
        p1.Rapidity(), p2.Rapidity(), p3.Rapidity(), p4.Rapidity(),
        // Four Pt
        p1234.Pt(),
        // Four Eta
        p1234.Eta(),
        // Four Phi
        p1234.Phi(),
        // Four Rapidity
        p1234.Rapidity(),
        // Four Mass
        p1234.M(),
        // Four Collins Soper Phi and CosTheta
        fourPiPhiPair1, fourPiPhiPair2, fourPiCosThetaPair1, fourPiCosThetaPair2);

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_0_charge_within_rap"), p1234.Pt());
        histosData.fill(HIST("fourpion_eta_0_charge_within_rap"), p1234.Eta());
        histosData.fill(HIST("fourpion_phi_0_charge_within_rap"), p1234.Phi());
        histosData.fill(HIST("fourpion_rap_0_charge_within_rap"), p1234.Rapidity());
        histosData.fill(HIST("fourpion_mass_0_charge_within_rap"), p1234.M());
        if (p1234.Pt() < rhoPtCut) {
          histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 14);
          if (rhoMassMin < p1234.M() && p1234.M() < rhoMassMax) {
            histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 15);
          }
          // Fill the Invariant Mass Histogram
          histosData.fill(HIST("fourpion_mass_0_charge_domA"), p1234.M());
          // Two Pion Masses
          histosData.fill(HIST("twopion_mass_1"), (p1 + p3).M());
          histosData.fill(HIST("twopion_mass_2"), (p1 + p4).M());
          histosData.fill(HIST("twopion_mass_3"), (p2 + p3).M());
          histosData.fill(HIST("twopion_mass_4"), (p2 + p4).M());
          // Fill the Collins-Soper Frame histograms
          histosData.fill(HIST("collin_soper_phi_1"), fourPiPhiPair1);
          histosData.fill(HIST("collin_soper_phi_2"), fourPiPhiPair2);
          histosData.fill(HIST("collin_soper_costheta_1"), fourPiCosThetaPair1);
          histosData.fill(HIST("collin_soper_costheta_2"), fourPiCosThetaPair2);
          histosData.fill(HIST("phi_vs_costheta_1"), fourPiPhiPair1, fourPiCosThetaPair1);
          histosData.fill(HIST("phi_vs_costheta_2"), fourPiPhiPair2, fourPiCosThetaPair2);

          // Small Mass CosTheta and Phi
          if ((k13.M() + k24.M()) > (k14.M() + k23.M())) {
            histosData.fill(HIST("collin_soper_phi_large_mass"), fourPiPhiPair1);
            histosData.fill(HIST("collin_soper_costheta_large_mass"), fourPiCosThetaPair1);
            histosData.fill(HIST("phi_vs_costheta_large_mass"), fourPiPhiPair1, fourPiCosThetaPair1);
            histosData.fill(HIST("collin_soper_phi_small_mass"), fourPiPhiPair2);
            histosData.fill(HIST("collin_soper_costheta_small_mass"), fourPiCosThetaPair2);
            histosData.fill(HIST("phi_vs_costheta_small_mass"), fourPiPhiPair2, fourPiCosThetaPair2);
          } else {
            histosData.fill(HIST("collin_soper_phi_small_mass"), fourPiPhiPair1);
            histosData.fill(HIST("collin_soper_costheta_small_mass"), fourPiCosThetaPair1);
            histosData.fill(HIST("phi_vs_costheta_small_mass"), fourPiPhiPair1, fourPiCosThetaPair1);
            histosData.fill(HIST("collin_soper_phi_large_mass"), fourPiPhiPair2);
            histosData.fill(HIST("collin_soper_costheta_large_mass"), fourPiCosThetaPair2);
            histosData.fill(HIST("phi_vs_costheta_large_mass"), fourPiPhiPair2, fourPiCosThetaPair2);
          }
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosData.fill(HIST("fourpion_mass_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosData.fill(HIST("fourpion_mass_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection
    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTracks != numPiMinus && numPiPlusTracks != numPiPlus) {

      histosData.fill(HIST("EventsCounts_vs_runNo"), collision.runNumber(), 16);

      ROOT::Math::PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p1234 = p1 + p2 + p3 + p4;

      histosData.fill(HIST("fourpion_pT_non_0_charge"), p1234.Pt());
      histosData.fill(HIST("fourpion_eta_non_0_charge"), p1234.Eta());
      histosData.fill(HIST("fourpion_phi_non_0_charge"), p1234.Phi());
      histosData.fill(HIST("fourpion_rap_non_0_charge"), p1234.Rapidity());
      histosData.fill(HIST("fourpion_mass_non_0_charge"), p1234.M());

      bkgFromData(
        // Run Number
        collision.runNumber(),
        // UPC mode
        collision.flags(),
        // vertex
        collision.posX(), collision.posY(), collision.posZ(),
        // FIT Signals
        collision.totalFV0AmplitudeA(), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
        // FIT and ZDC Signals
        collision.timeFV0A(), collision.timeFT0A(), collision.timeFT0C(), collision.timeFDDA(), collision.timeFDDC(), collision.timeZNA(), collision.timeZNC(),
        // Occupancy
        collision.occupancyInTime(),
        // has atleast one TOF
        hasAtleastOneTOF,
        // DCA XY and Z
        selectedPionTracks[0].dcaXY(), selectedPionTracks[1].dcaXY(), selectedPionTracks[2].dcaXY(), selectedPionTracks[3].dcaXY(),
        selectedPionTracks[0].dcaZ(), selectedPionTracks[1].dcaZ(), selectedPionTracks[2].dcaZ(), selectedPionTracks[3].dcaZ(),
        // TPC N Sigma Pi
        selectedPionTracks[0].tpcNSigmaPi(), selectedPionTracks[1].tpcNSigmaPi(), selectedPionTracks[2].tpcNSigmaPi(), selectedPionTracks[3].tpcNSigmaPi(),
        // TPC N Sigma Ka
        selectedPionTracks[0].tpcNSigmaKa(), selectedPionTracks[1].tpcNSigmaKa(), selectedPionTracks[2].tpcNSigmaKa(), selectedPionTracks[3].tpcNSigmaKa(),
        // TPC N Sigma Pr
        selectedPionTracks[0].tpcNSigmaPr(), selectedPionTracks[1].tpcNSigmaPr(), selectedPionTracks[2].tpcNSigmaPr(), selectedPionTracks[3].tpcNSigmaPr(),
        // TPC N Sigma El
        selectedPionTracks[0].tpcNSigmaEl(), selectedPionTracks[1].tpcNSigmaEl(), selectedPionTracks[2].tpcNSigmaEl(), selectedPionTracks[3].tpcNSigmaEl(),
        // TPC N Sigma Mu
        selectedPionTracks[0].tpcNSigmaMu(), selectedPionTracks[1].tpcNSigmaMu(), selectedPionTracks[2].tpcNSigmaMu(), selectedPionTracks[3].tpcNSigmaMu(),
        // tpc Chi2 NCl
        selectedPionTracks[0].tpcChi2NCl(), selectedPionTracks[1].tpcChi2NCl(), selectedPionTracks[2].tpcChi2NCl(), selectedPionTracks[3].tpcChi2NCl(),
        // TPC NCls Findable
        selectedPionTracks[0].tpcNClsFindable(), selectedPionTracks[1].tpcNClsFindable(), selectedPionTracks[2].tpcNClsFindable(), selectedPionTracks[3].tpcNClsFindable(),
        // ITS Chi2 NCl
        selectedPionTracks[0].itsChi2NCl(), selectedPionTracks[1].itsChi2NCl(), selectedPionTracks[2].itsChi2NCl(), selectedPionTracks[3].itsChi2NCl(),
        // Pion Pt
        p1.Pt(), p2.Pt(), p3.Pt(), p4.Pt(),
        // Pion Eta
        p1.Eta(), p2.Eta(), p3.Eta(), p4.Eta(),
        // Pion Phi
        p1.Phi(), p2.Phi(), p3.Phi(), p4.Phi(),
        // Pion Rapidity
        p1.Rapidity(), p2.Rapidity(), p3.Rapidity(), p4.Rapidity(),
        // Four Pt
        p1234.Pt(),
        // Four Eta
        p1234.Eta(),
        // Four Phi
        p1234.Phi(),
        // Four Rapidity
        p1234.Rapidity(),
        // Four Mass
        p1234.M());

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_non_0_charge_within_rap"), p1234.Pt());
        histosData.fill(HIST("fourpion_eta_non_0_charge_within_rap"), p1234.Eta());
        histosData.fill(HIST("fourpion_phi_non_0_charge_within_rap"), p1234.Phi());
        histosData.fill(HIST("fourpion_rap_non_0_charge_within_rap"), p1234.Rapidity());
        histosData.fill(HIST("fourpion_mass_non_0_charge_within_rap"), p1234.M());
        if (p1234.Pt() < rhoPtCut) {
          histosData.fill(HIST("fourpion_mass_non_0_charge_domA"), p1234.M());
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosData.fill(HIST("fourpion_mass_non_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosData.fill(HIST("fourpion_mass_non_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection
    } // End of Analysis for non 0 charge events
  } // End of 4 Pion Analysis Process function for Pass5 Data

  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processData, "The Process for 4 Pion Analysis from data", true);

  double cosThetaCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);
    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;
    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2CM{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};
    double cosThetaCS = zaxisCS.Dot((v1CM));
    return cosThetaCS;
  }

  double phiCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);
    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;
    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2CM{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCS{((beam1CM.Unit() - beam2CM.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxisCS{(beam1CM.Cross(beam2CM)).Unit()};
    ROOT::Math::XYZVectorF xaxisCS{(yaxisCS.Cross(zaxisCS)).Unit()};

    double phi = std::atan2(yaxisCS.Dot(v1CM), xaxisCS.Dot(v1CM));
    return phi;
  }

  double getMaxDCAxy(double pT)
  {
    return 0.0105 + 0.035 / std::pow(pT, 1.1);
  }

}; // End of Struct exclusiveRhoTo4Pi

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusiveRhoTo4Pi>(cfgc)};
}
