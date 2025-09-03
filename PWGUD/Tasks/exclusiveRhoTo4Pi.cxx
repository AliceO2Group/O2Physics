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
// Collin-Soper Angles
DECLARE_SOA_COLUMN(FourPionPhiPair1, fourPionPhiPair1, double);
DECLARE_SOA_COLUMN(FourPionPhiPair2, fourPionPhiPair2, double);
DECLARE_SOA_COLUMN(FourPionPhiPair3, fourPionPhiPair3, double);
DECLARE_SOA_COLUMN(FourPionPhiPair4, fourPionPhiPair4, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair1, fourPionCosThetaPair1, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair2, fourPionCosThetaPair2, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair3, fourPionCosThetaPair3, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair4, fourPionCosThetaPair4, double);
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
                  branch::FourPionPhiPair3,
                  branch::FourPionPhiPair4,
                  branch::FourPionCosThetaPair1,
                  branch::FourPionCosThetaPair2,
                  branch::FourPionCosThetaPair3,
                  branch::FourPionCosThetaPair4);

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
  double mRho0 = 0.77526; // GeV/c^2
  // Run Numbers
  static int runNos[113];
  static int numRunNums;
  // Derived Data
  Produces<aod::SignalData> sigFromData;
  Produces<aod::BkgroundData> bkgFromData;
  // Histogram Registry
  HistogramRegistry histosData{"Data", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosCounter{"counters", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  // Configurable Event parameters
  Configurable<int> ifUPC{"ifUPC", 1, "Enable UPC reconstruction only"};
  Configurable<float> vZCut{"vZCut", 10., "Vertex Cut"};
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 50., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> zdcCut{"zdcCut", 0., "ZDC threshold"};
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
  Configurable<int> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min TPC Findable Clusters"};
  // Configurable PID parameters
  Configurable<bool> useTOF{"useTOF", true, "if track has TOF use TOF"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 3, "TPC cut"};
  Configurable<float> nSigmaTOFcut{"nSigmaTOFcut", 3, "TOF cut"};
  // Configurable Rho parameters
  Configurable<float> rhoRapCut{"rhoRapCut", 0.5, "Max abs Rapidity of rho"};
  Configurable<float> rhoPtCut{"rhoPtCut", 0.15, "Min Pt of rho"};
  Configurable<float> rhoMassMin{"rhoMassMin", 1, "Min Mass of rho"};
  Configurable<float> rhoMassMax{"rhoMassMax", 2.5, "Max Mass of rho"};
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
    histosCounter.add("EventsCounts_vs_runNo", "Event Counter Run by Run; Run Number; Number of Events", kTH2F, {{113, 0, 113}, {14, 0, 14}});
    histosCounter.add("TracksCounts_vs_runNo", "Track Counter Run by Run; Run Number; Number of Tracks", kTH2F, {{113, 0, 113}, {14, 0, 14}});
    histosCounter.add("fourPionCounts_0c", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_0c_within_mass", "Four Pion Counts within mass range; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_0c_within_rap", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_0c_selected", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_n0c", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_n0c_within_mass", "Four Pion Counts within mass range; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_n0c_within_rap", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    histosCounter.add("fourPionCounts_n0c_selected", "Four Pion Counts; Run Number; Events", kTH1F, {{113, 0, 113}});
    // QA plots: event selection
    histosData.add("UPCmode", "UPC mode; Events", kTH1F, {{5, 0, 5}});
    histosData.add("GapSide", "Gap Side;Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosData.add("TrueGapSide", "True Gap Side; True Gap Side; Events", kTH1F, {{4, 0, 4}});
    histosData.add("isCBTOk", "isCBTOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosData.add("isCBTHadronOk", "isCBTHadronOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosData.add("isCBTZdcOk", "isCBTZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosData.add("isCBTHadronZdcOk", "isCBTHadronZdcOk; bool; Events", kTH1F, {{4, 0, 4}});
    histosData.add("FT0A", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosData.add("FT0C", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
    histosData.add("FV0A", "V0A amplitude", kTH1F, {{100, 0.0, 100}});
    histosData.add("ZDC_A", "ZDC amplitude", kTH1F, {{10000, 0.0, 10000}});
    histosData.add("ZDC_C", "ZDC amplitude", kTH1F, {{10000, 0.0, 10000}});
    histosData.add("FDDA", "FDD A signal; FDD A signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosData.add("FDDC", "FDD C signal; FDD C signal; Counts", kTH1F, {{500, 0.0, 2000}});
    histosData.add("vertexX", "Vertex X; Vertex X [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosData.add("vertexY", "Vertex Y; Vertex Y [cm]; Counts", kTH1F, {{2000, -0.05, 0.05}});
    histosData.add("vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{2000, -15, 15}});
    histosData.add("occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{20000, 0, 20000}});
    // QA plots: tracks
    histosData.add("dcaXY_all", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosData.add("dcaXY_pions", "dcaXY_pions; dcaXY of Pions [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosData.add("dcaZ_all", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosData.add("dcaZ_pions", "dcaZ_pions; dcaZ of Pions [cm]; Counts", kTH1F, {{2000, -0.1, 0.1}});
    histosData.add("itsChi2NCl_all", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("itsChi2_all", "ITS Chi2; ITS Chi2; Counts", kTH1F, {{500, 0, 50}});
    histosData.add("tpcChi2NCl_all", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("tpcNClsFindable_all", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});
    // QA plots: PID
    histosData.add("tpcSignal_all", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosData.add("tpcSignal_pions", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosData.add("tpcNSigmaPi_all", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPi_pions", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaKa_pions", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPr_pions", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaEl_pions", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaMu_pions", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofBeta_all", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
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
    // Collin Soper Theta and Phi after selection
    histosData.add("CSphi_vs_CScosTheta", "Phi vs cosTheta for small mass; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    setHistBinLabels();
  } // End of init function

  //---------------------------------------------------------------------------------------------------------------------------------------------
  // Event Cuts
  Filter vertexZcut = (nabs(o2::aod::collision::posZ) <= vZCut);
  Filter numPVcontributorsCut = (o2::aod::collision::numContrib == numPVContrib);
  Filter fitcuts = (o2::aod::udcollision::totalFV0AmplitudeA <= fv0Cut) && (o2::aod::udcollision::totalFT0AmplitudeA <= ft0aCut) && (o2::aod::udcollision::totalFT0AmplitudeC <= ft0cCut);
  Filter zdcCuts = (o2::aod::udzdc::energyCommonZNA <= zdcCut) && (o2::aod::udzdc::energyCommonZNC <= zdcCut);
  Filter bcSelectionCuts = (o2::aod::udcollision::sbp == sbpCut) && (o2::aod::udcollision::itsROFb == itsROFbCut) && (o2::aod::udcollision::vtxITSTPC == vtxITSTPCcut) && (o2::aod::udcollision::tfb == tfbCut);
  // Track Cuts
  Filter onlyPVtracks = o2::aod::udtrack::isPVContributor == useOnlyPVtracks;
  Filter tpcchi2nclsFilter = o2::aod::track::tpcChi2NCl <= tpcChi2NClsCut;
  Filter itschi2nclsFilter = o2::aod::track::itsChi2NCl <= itsChi2NClsCut;
  Filter tpcCuts = (nabs(o2::aod::pidtpc::tpcNSigmaPi) <= nSigmaTPCcut);
  //---------------------------------------------------------------------------------------------------------------------------------------------

  using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;

  void processData(soa::Filtered<UDCollisions>::iterator const& collision, soa::Filtered<UDtracks> const& tracks)
  {

    // Check if the Event is reconstructed in UPC mode
    if (collision.flags() != ifUPC) {
      return;
    }

    // RCT flag
    if (!sgSelector.isCBTHadronZdcOk(collision)) {
      return;
    }

    int runIndex = getRunNumberIndex(collision.runNumber());

    histosData.fill(HIST("GapSide"), collision.gapSide());
    histosData.fill(HIST("TrueGapSide"), sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut));
    histosData.fill(HIST("isCBTOk"), sgSelector.isCBTOk(collision));
    histosData.fill(HIST("isCBTHadronOk"), sgSelector.isCBTHadronOk(collision));
    histosData.fill(HIST("isCBTZdcOk"), sgSelector.isCBTZdcOk(collision));
    histosData.fill(HIST("isCBTHadronZdcOk"), sgSelector.isCBTHadronZdcOk(collision));
    histosData.fill(HIST("vertexX"), collision.posX());
    histosData.fill(HIST("vertexY"), collision.posY());
    histosData.fill(HIST("vertexZ"), collision.posZ());
    histosData.fill(HIST("occupancy"), collision.occupancyInTime());
    histosData.fill(HIST("FV0A"), collision.totalFV0AmplitudeA());
    histosData.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histosData.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histosData.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histosData.fill(HIST("ZDC_C"), collision.energyCommonZNC());
    histosData.fill(HIST("FDDA"), collision.totalFDDAmplitudeA());
    histosData.fill(HIST("FDDC"), collision.totalFDDAmplitudeC());
    histosData.fill(HIST("UPCmode"), collision.flags());

    std::vector<decltype(tracks.begin())> selectedTracks;
    std::vector<decltype(tracks.begin())> selectedPionTracks;
    std::vector<decltype(tracks.begin())> selectedPionPlusTracks;
    std::vector<decltype(tracks.begin())> selectedPionMinusTracks;

    for (const auto& t0 : tracks) {
      if (!isSelectedTrack(t0, pTcut, etaCut, dcaXYcut, dcaZcut, useITStracksOnly, useTPCtracksOnly, itsChi2NClsCut, tpcChi2NClsCut, tpcNClsFindableCut)) {
        continue;
      }
      selectedTracks.push_back(t0);
      if (selectionPIDPion(t0, useTOF, nSigmaTPCcut, nSigmaTOFcut)) {
        selectedPionTracks.push_back(t0);
        if (t0.sign() == 1) {
          selectedPionPlusTracks.push_back(t0);
        }
        if (t0.sign() == -1) {
          selectedPionMinusTracks.push_back(t0);
        }
      } // End of Selection PID Pion
    } // End of loop over tracks

    int numSelectedTracks = static_cast<int>(selectedTracks.size());
    int numSelectedPionTracks = static_cast<int>(selectedPionTracks.size());
    int numPiPlusTracks = static_cast<int>(selectedPionPlusTracks.size());
    int numPionMinusTracks = static_cast<int>(selectedPionMinusTracks.size());

    for (int i = 0; i < numSelectedTracks; i++) {
      PxPyPzMVector selectedTrackVector(selectedTracks[i].px(), selectedTracks[i].py(), selectedTracks[i].pz(), o2::constants::physics::MassPionCharged);
      histosData.fill(HIST("pT_track_all"), selectedTrackVector.Pt());
      histosData.fill(HIST("eta_track_all"), selectedTrackVector.Eta());
      histosData.fill(HIST("phi_track_all"), selectedTrackVector.Phi());
      histosData.fill(HIST("rapidity_track_all"), selectedTrackVector.Rapidity());

      histosData.fill(HIST("dcaXY_all"), selectedTracks[i].dcaXY());
      histosData.fill(HIST("dcaZ_all"), selectedTracks[i].dcaZ());

      histosData.fill(HIST("itsChi2NCl_all"), selectedTracks[i].itsChi2NCl());
      histosData.fill(HIST("itsChi2_all"), selectedTracks[i].itsChi2NCl() * selectedTracks[i].itsNCls());
      histosData.fill(HIST("tpcChi2NCl_all"), selectedTracks[i].tpcChi2NCl());
      histosData.fill(HIST("tpcNClsFindable_all"), selectedTracks[i].tpcNClsFindable());

      histosData.fill(HIST("tpcSignal_all"), selectedTrackVector.P(), selectedTracks[i].tpcSignal());
      histosData.fill(HIST("tpcNSigmaPi_all"), selectedTracks[i].tpcNSigmaPi(), selectedTrackVector.Pt());
      histosData.fill(HIST("tofBeta_all"), selectedTrackVector.P(), selectedTracks[i].beta());
      histosData.fill(HIST("tofNSigmaPi_all"), selectedTracks[i].tofNSigmaPi(), selectedTrackVector.Pt());
    } // End of loop over tracks with selection only

    for (int i = 0; i < numSelectedPionTracks; i++) {
      PxPyPzMVector selectedPionTrackVector(selectedPionTracks[i].px(), selectedPionTracks[i].py(), selectedPionTracks[i].pz(), o2::constants::physics::MassPionCharged);

      histosData.fill(HIST("pT_track_pions"), selectedPionTrackVector.Pt());
      histosData.fill(HIST("eta_track_pions"), selectedPionTrackVector.Eta());
      histosData.fill(HIST("phi_track_pions"), selectedPionTrackVector.Phi());
      histosData.fill(HIST("rapidity_track_pions"), selectedPionTrackVector.Rapidity());

      histosData.fill(HIST("dcaXY_pions"), selectedPionTracks[i].dcaXY());
      histosData.fill(HIST("dcaZ_pions"), selectedPionTracks[i].dcaZ());

      histosData.fill(HIST("tpcSignal_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].tpcSignal());
      histosData.fill(HIST("tpcNSigmaPi_pions"), selectedPionTracks[i].tpcNSigmaPi(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaKa_pions"), selectedPionTracks[i].tpcNSigmaKa(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaPr_pions"), selectedPionTracks[i].tpcNSigmaPr(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaEl_pions"), selectedPionTracks[i].tpcNSigmaEl(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tpcNSigmaMu_pions"), selectedPionTracks[i].tpcNSigmaMu(), selectedPionTrackVector.Pt());

      histosData.fill(HIST("tofBeta_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].beta());
      histosData.fill(HIST("tofNSigmaPi_pions"), selectedPionTracks[i].tofNSigmaPi(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaKa_pions"), selectedPionTracks[i].tofNSigmaKa(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaPr_pions"), selectedPionTracks[i].tofNSigmaPr(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaEl_pions"), selectedPionTracks[i].tofNSigmaEl(), selectedPionTrackVector.Pt());
      histosData.fill(HIST("tofNSigmaMu_pions"), selectedPionTracks[i].tofNSigmaMu(), selectedPionTrackVector.Pt());
    } // End of loop over tracks with selection and PID of pions

    // event should have exactly 4 pions
    if (numSelectedPionTracks != numFourPionTracks) {
      return;
    }

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

      PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

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

      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;
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

      double fourPiPhiPair1 = collinSoperPhi(k13, k1234);
      double fourPiPhiPair2 = collinSoperPhi(k14, k1234);
      double fourPiPhiPair3 = collinSoperPhi(k23, k1234);
      double fourPiPhiPair4 = collinSoperPhi(k24, k1234);
      double fourPiCosThetaPair1 = collinSoperCosTheta(k13, k1234);
      double fourPiCosThetaPair2 = collinSoperCosTheta(k14, k1234);
      double fourPiCosThetaPair3 = collinSoperCosTheta(k23, k1234);
      double fourPiCosThetaPair4 = collinSoperCosTheta(k24, k1234);

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
        fourPiPhiPair1, fourPiPhiPair2, fourPiPhiPair3, fourPiPhiPair4,
        fourPiCosThetaPair1, fourPiCosThetaPair2, fourPiCosThetaPair3, fourPiCosThetaPair4);

      histosCounter.fill(HIST("fourPionCounts_0c"), runIndex);

      if (rhoMassMin < p1234.M() && p1234.M() < rhoMassMax) {
        histosCounter.fill(HIST("fourPionCounts_0c_within_mass"), runIndex);
      }

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_0_charge_within_rap"), p1234.Pt());
        histosData.fill(HIST("fourpion_eta_0_charge_within_rap"), p1234.Eta());
        histosData.fill(HIST("fourpion_phi_0_charge_within_rap"), p1234.Phi());
        histosData.fill(HIST("fourpion_rap_0_charge_within_rap"), p1234.Rapidity());
        histosData.fill(HIST("fourpion_mass_0_charge_within_rap"), p1234.M());
        histosCounter.fill(HIST("fourPionCounts_0c_within_rap"), runIndex);
        if (p1234.Pt() < rhoPtCut) {
          if (rhoMassMin < p1234.M() && p1234.M() < rhoMassMax) {
            // Selected Four Pion Events
            histosCounter.fill(HIST("fourPionCounts_0c_selected"), runIndex);
            // Fill the Invariant Mass Histogram
            histosData.fill(HIST("fourpion_mass_0_charge_domA"), p1234.M());
            // Two Pion Masses
            histosData.fill(HIST("twopion_mass_1"), (p1 + p3).M());
            histosData.fill(HIST("twopion_mass_2"), (p1 + p4).M());
            histosData.fill(HIST("twopion_mass_3"), (p2 + p3).M());
            histosData.fill(HIST("twopion_mass_4"), (p2 + p4).M());
            // Fill the Collins-Soper Frame histograms
            double mDiff13 = std::abs((k13.M() - mRho0));
            double mDiff14 = std::abs((k14.M() - mRho0));
            double mDiff23 = std::abs((k23.M() - mRho0));
            double mDiff24 = std::abs((k24.M() - mRho0));
            if ((mDiff13 < mDiff14) && (mDiff13 < mDiff23) && (mDiff13 < mDiff24)) {
              histosData.fill(HIST("CSphi_vs_CScosTheta"), fourPiPhiPair1, fourPiCosThetaPair1);
            } else if ((mDiff14 < mDiff13) && (mDiff14 < mDiff23) && (mDiff14 < mDiff24)) {
              histosData.fill(HIST("CSphi_vs_CScosTheta"), fourPiPhiPair2, fourPiCosThetaPair2);
            } else if ((mDiff23 < mDiff13) && (mDiff23 < mDiff14) && (mDiff23 < mDiff24)) {
              histosData.fill(HIST("CSphi_vs_CScosTheta"), fourPiPhiPair3, fourPiCosThetaPair3);
            } else if ((mDiff24 < mDiff13) && (mDiff24 < mDiff14) && (mDiff24 < mDiff23)) {
              histosData.fill(HIST("CSphi_vs_CScosTheta"), fourPiPhiPair4, fourPiCosThetaPair4);
            }
          } // End of Pt selection for rho mass
        } // End of Pt selection for rho mass
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

      PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);
      PxPyPzMVector p1234 = p1 + p2 + p3 + p4;

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

      histosCounter.fill(HIST("fourPionCounts_n0c"), runIndex);
      if (rhoMassMin < p1234.M() && p1234.M() < rhoMassMax) {
        histosCounter.fill(HIST("fourPionCounts_n0c_within_mass"), runIndex);
      }
      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_non_0_charge_within_rap"), p1234.Pt());
        histosData.fill(HIST("fourpion_eta_non_0_charge_within_rap"), p1234.Eta());
        histosData.fill(HIST("fourpion_phi_non_0_charge_within_rap"), p1234.Phi());
        histosData.fill(HIST("fourpion_rap_non_0_charge_within_rap"), p1234.Rapidity());
        histosData.fill(HIST("fourpion_mass_non_0_charge_within_rap"), p1234.M());
        histosCounter.fill(HIST("fourPionCounts_n0c_within_rap"), runIndex);
        if (p1234.Pt() < rhoPtCut) {
          histosData.fill(HIST("fourpion_mass_non_0_charge_domA"), p1234.M());
          histosCounter.fill(HIST("fourPionCounts_n0c_selected"), runIndex);
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

  void processEventCounter(UDCollisions::iterator const& collision)
  {
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 0);
    // RCT flag
    if (!sgSelector.isCBTHadronZdcOk(collision)) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 1);
    // UPC mode
    if (collision.flags() != ifUPC) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 2);
    // vtxITSTPC
    if (collision.vtxITSTPC() != vtxITSTPCcut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 3);
    // sbp
    if (collision.sbp() != sbpCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 4);
    // itsROFb
    if (collision.itsROFb() != itsROFbCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 5);
    // tfb
    if (collision.tfb() != tfbCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 6);
    // FT0A
    if (collision.totalFT0AmplitudeA() > ft0aCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 7);
    // FT0C
    if (collision.totalFT0AmplitudeC() > ft0cCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 8);
    // FV0A
    if (collision.totalFV0AmplitudeA() > fv0Cut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 9);
    // ZDC
    if (collision.energyCommonZNA() > zdcCut || collision.energyCommonZNC() > zdcCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 10);
    // numContributors
    if (collision.numContrib() != numPVContrib) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 11);
    // vertexZ
    if (std::abs(collision.posZ()) > vZCut) {
      return;
    }
    histosCounter.fill(HIST("EventsCounts_vs_runNo"), getRunNumberIndex(collision.runNumber()), 12);
  } // End of processCounter function

  void processTrackCounter(soa::Filtered<UDCollisions>::iterator const& collision, UDtracks const& tracks)
  {
    int runIndex = getRunNumberIndex(collision.runNumber());
    // Check if the Event is reconstructed in UPC mode
    if (collision.flags() != ifUPC) {
      return;
    }
    for (const auto& track : tracks) {
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 0);
      PxPyPzMVector trackVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);
      // is PV contributor
      if (track.isPVContributor() != useOnlyPVtracks) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 1);
      // pt cut
      if (trackVector.Pt() < pTcut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 2);
      // eta cut
      if (std::abs(trackVector.Eta()) > etaCut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 3);
      // DCA Z cut
      if (std::abs(track.dcaZ()) > dcaZcut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 4);
      // DCA XY cut
      float maxDCAxy = 0.0105 + 0.035 / std::pow(trackVector.Pt(), 1.1);
      if (dcaXYcut == 0 && (std::fabs(track.dcaXY()) > maxDCAxy)) {
        continue;
      } else if (dcaXYcut != 0 && (std::fabs(track.dcaXY()) > dcaXYcut)) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 5);
      // ITS Track only
      if (useITStracksOnly && !track.hasITS()) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 6);
      // TPC Track only
      if (useTPCtracksOnly && !track.hasTPC()) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 7);
      // ITS Chi2 N Clusters cut
      if (track.hasITS() && track.itsChi2NCl() > itsChi2NClsCut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 8);
      // TPC Chi2 N Clusters cut
      if (track.hasTPC() && track.tpcChi2NCl() > tpcChi2NClsCut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 9);
      // TPC N Clusters Findable cut
      if (track.hasTPC() && track.tpcNClsFindable() < tpcNClsFindableCut) {
        continue;
      }
      histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 10);
      // Selection PID Pion
      if (selectionPIDPion(track, useTOF, nSigmaTPCcut, nSigmaTOFcut)) {
        histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 11);
        if (track.sign() == 1) {
          histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 12);
        }
        if (track.sign() == -1) {
          histosCounter.fill(HIST("TracksCounts_vs_runNo"), runIndex, 13);
        }
      } // End of Selection PID Pion
    } // End of loop over tracks
  } // End of processCounter function

  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processData, "Data Analysis Function", true);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processEventCounter, "Event Counter Function", true);
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processTrackCounter, "Track Counter Function", true);

  double collinSoperPhi(PtEtaPhiMVector twoPionVector, PtEtaPhiMVector fourPionVector)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);
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

  double collinSoperCosTheta(PtEtaPhiMVector twoPionVector, PtEtaPhiMVector fourPionVector)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);
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
                       float tpcnclsfindablecut)
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
    // ITS Chi2 N Clusters cut
    if (track.hasITS() && track.itsChi2NCl() > itschi2nclscut) {
      return false;
    }
    // TPC Chi2 N Clusters cut
    if (track.hasTPC() && track.tpcChi2NCl() > tpcchi2nclscut) {
      return false;
    }
    // TPC N Clusters Findable cut
    if (track.hasTPC() && track.tpcNClsFindable() < tpcnclsfindablecut) {
      return false;
    }
    // All cuts passed
    return true;
  } // End of Track Selection function

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
      "ZDC",
      "n PV Contrib = 4",
      "V_{z} < " + strFormat(vZCut) + " cm"};

    int numEventCuts = 13;

    std::string trackLabels[14] = {
      "No Cuts",
      "isPVContributor",
      "pT>" + strFormat(pTcut) + " GeV/c",
      "|#eta|<" + strFormat(etaCut),
      "DCA Z<" + strFormat(dcaZcut) + " cm",
      "DCA XY cut",
      "hasITS",
      "hasTPC",
      "itsChi2NCl<" + strFormat(itsChi2NClsCut),
      "tpcChi2NCl<" + strFormat(tpcChi2NClsCut),
      "tpcNClsFindable>" + strFormat(tpcNClsFindableCut),
      "#pi tracks",
      "#pi^{+} tracks",
      "#pi^{-} tracks"};

    int numTrackCuts = 14;

    auto h1 = histosCounter.get<TH2>(HIST("EventsCounts_vs_runNo"));
    auto h2 = histosCounter.get<TH2>(HIST("TracksCounts_vs_runNo"));
    auto h3 = histosCounter.get<TH1>(HIST("fourPionCounts_0c"));
    auto h4 = histosCounter.get<TH1>(HIST("fourPionCounts_0c_within_rap"));
    auto h5 = histosCounter.get<TH1>(HIST("fourPionCounts_0c_selected"));
    auto h6 = histosCounter.get<TH1>(HIST("fourPionCounts_n0c"));
    auto h7 = histosCounter.get<TH1>(HIST("fourPionCounts_n0c_within_rap"));
    auto h8 = histosCounter.get<TH1>(HIST("fourPionCounts_n0c_selected"));
    auto h9 = histosCounter.get<TH1>(HIST("fourPionCounts_0c_within_mass"));
    auto h10 = histosCounter.get<TH1>(HIST("fourPionCounts_n0c_within_mass"));

    for (int i = 0; i < numRunNums; ++i) {
      h1->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h2->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h3->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h4->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h5->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h6->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h7->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h8->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h9->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
      h10->GetXaxis()->SetBinLabel(i + 1, std::to_string(runNos[i]).c_str());
    }
    for (int i = 0; i < numEventCuts; ++i) {
      h1->GetYaxis()->SetBinLabel(i + 1, eventLabels[i].c_str());
    }
    for (int i = 0; i < numTrackCuts; ++i) {
      h2->GetYaxis()->SetBinLabel(i + 1, trackLabels[i].c_str());
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
