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

// DCA
DECLARE_SOA_COLUMN(Dcaxy1, dcaxy1, double);
DECLARE_SOA_COLUMN(Dcaxy2, dcaxy2, double);
DECLARE_SOA_COLUMN(Dcaxy3, dcaxy3, double);
DECLARE_SOA_COLUMN(Dcaxy4, dcaxy4, double);

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

DECLARE_SOA_COLUMN(FourPionPt, fourPionPt, double);
DECLARE_SOA_COLUMN(FourPionEta, fourPionEta, double);
DECLARE_SOA_COLUMN(FourPionPhi, fourPionPhi, double);
DECLARE_SOA_COLUMN(FourPionRapidity, fourPionRapidity, double);

DECLARE_SOA_COLUMN(FourPionMass, fourPionMass, double);
DECLARE_SOA_COLUMN(FourPionPhiPair1, fourPionPhiPair1, double);
DECLARE_SOA_COLUMN(FourPionPhiPair2, fourPionPhiPair2, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair1, fourPionCosThetaPair1, double);
DECLARE_SOA_COLUMN(FourPionCosThetaPair2, fourPionCosThetaPair2, double);
} // namespace branch

DECLARE_SOA_TABLE(SignalData, "AOD", "signalData",
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

DECLARE_SOA_TABLE(MCgen, "AOD", "MCgen",
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

DECLARE_SOA_TABLE(SignalMCreco, "AOD", "SignalMCreco",
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
} // namespace o2::aod

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct ExclusiveRhoTo4Pi {
  SGSelector sgSelector;
  int rhoPrime = 30113;
  int numFourPionTracks = 4;
  int numPiPlus = 2;
  int numPiMinus = 2;
  float zeroPointEight = 0.8;
  Produces<aod::SignalData> sigFromData;
  Produces<aod::BkgroundData> bkgFromData;
  Produces<aod::MCgen> generatedMC;
  Produces<aod::SignalMCreco> sigFromMC;

  HistogramRegistry histosData{"histosData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosMCgen{"histosMCgen", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosMCreco{"histosMCreco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Event Selection Cuts
  Configurable<float> vZCut{"vZCut", 10., "Vertex Cut"};
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 150., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> zdcCut{"zdcCut", 0., "ZDC threshold"};
  Configurable<float> occupancyCut{"occupancyCut", 20000, "Occupancy Cut"};
  Configurable<uint16_t> numPVContrib{"numPVContrib", 4, "Number of PV Contributors"};
  Configurable<int> checkOneTof{"checkOneTof", 1, " (1 or 0)Check if event has at least 1 TOF"};

  // bc selection cuts
  Configurable<int> sbpCut{"sbpCut", 1, "Sbp"};
  Configurable<int> itsROFbCut{"itsROFbCut", 1, "itsROFbCut"};
  Configurable<int> vtxITSTPCcut{"vtxITSTPCcut", 1, "vtxITSTPCcut"};
  Configurable<int> tfbCut{"tfbCut", 1, "tfbCut"};

  // Track Selection Cuts
  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZcut{"dcaZcut", 2, "dcaZ cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0, "dcaXY cut"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pTcut{"pTcut", 0.15, "Track Pt"};

  // PID Selection Cuts
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 3, "TPC cut"};
  Configurable<float> nSigmaTOFcut{"nSigmaTOFcut", 3, "TOF cut"};

  // Kinematics Cuts
  Configurable<float> rhoRapCut{"rhoRapCut", 0.5, "Max abs Rapidity of rho"};
  Configurable<float> rhoPtCut{"rhoPtCut", 0.15, "Min Pt of rho"};

  // Axis Configurations
  ConfigurableAxis pTAxis{"pTAxis", {1000, 0, 2}, "Axis for pT histograms"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {1000, -2.5, 2.5}, "Axis for Rapidity histograms"};
  ConfigurableAxis invMassAxis{"invMassAxis", {1000, 1, 2.5}, "Axis for Phi histograms"};
  ConfigurableAxis phiAxis{"phiAxis", {360, -1 * o2::constants::math::PI, o2::constants::math::PI}, "Axis for Phi histograms"};
  ConfigurableAxis cosThetaAxis{"cosThetaAxis", {360, -1, 1}, "Axis for cos Theta histograms"};
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Init Function-----------------------------------------------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {

    histosData.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histosData.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histosData.add("EventCounts", "Total Events; Events", kTH1F, {{10, 0, 10}});
    histosData.add("vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{2000, -15, 15}});
    histosData.add("occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{20000, 0, 20000}});
    histosData.add("dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{5000, -0.5, 0.5}});
    histosData.add("dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{5000, -0.5, 0.5}});
    histosData.add("tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{250, 0, 50}});
    histosData.add("tpcNClsFindable", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});

    // TPC nSigma
    histosData.add("tpcNSigmaPi_all", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPi_pions", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TPC nSigma of other particles with selected pion tracks
    histosData.add("tpcNSigmaKa_pions", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaPr_pions", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaEl_pions", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tpcNSigmaMu_pions", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TOF nSigma
    histosData.add("tofNSigmaPi_all", "TOF nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaPi_pions", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TOF nSigma of other particles with selected pion tracks
    histosData.add("tofNSigmaKa_pions", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaPr_pions", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaEl_pions", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosData.add("tofNSigmaMu_pions", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // Track Transverse Momentum
    histosData.add("pT_track_all", "pT with track selection; pT [GeV/c]; Counts", kTH1F, {pTAxis});
    histosData.add("pT_track_pions", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosData.add("pT_track_pions_contributed", "pT with track selection and PID selection of Pi which are contributed to selected event; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Track Rapidity
    histosData.add("rapidity_track_all", "Rapidity with track selection; y; Counts", kTH1F, {rapidityAxis});
    histosData.add("rapidity_track_pions", "Rapidity with track selection and PID selection of Pi; y; Events", kTH1F, {rapidityAxis});
    histosData.add("rapidity_track_pions_contributed", "Rapidity with track selection and PID selection of Pi which are contributed to selected event; y; Events", kTH1F, {rapidityAxis});

    // Zero charge Event Transverse Momentum
    histosData.add("fourpion_pT_0_charge", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Non Zero charge Event Transverse Momentum
    histosData.add("fourpion_pT_non_0_charge", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Rapidity of 0 charge Events
    histosData.add("fourpion_rap_0_charge_domA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosData.add("fourpion_rap_0_charge_domB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosData.add("fourpion_rap_0_charge_domC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});

    // Rapidity of non 0 charge Events
    histosData.add("fourpion_rap_non_0_charge_domA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {rapidityAxis});
    histosData.add("fourpion_rap_non_0_charge_domB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c$; y; Events", kTH1F, {rapidityAxis});
    histosData.add("fourpion_rap_non_0_charge_domC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {rapidityAxis});

    // Pair Invariant Mass
    histosData.add("twopion_mass_1", "Invariant Mass Distribution of 2 pions 1 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_2", "Invariant Mass Distribution of 2 pions 2 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_3", "Invariant Mass Distribution of 2 pions 3 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosData.add("twopion_mass_4", "Invariant Mass Distribution of 2 pions 4 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});

    // Invariant Mass of 0 charge events
    histosData.add("fourpion_mass_0_charge_domA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // pT < 0.15GeV
    histosData.add("fourpion_mass_0_charge_domB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis}); // 0.15GeV < pT < 0.8GeV
    histosData.add("fourpion_mass_0_charge_domC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // 0.8GeV < pT

    // Invariant mass of non 0 charge events
    histosData.add("fourpion_mass_non_0_charge_domA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // pT < 0.15GeV
    histosData.add("fourpion_mass_non_0_charge_domB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis}); // 0.15GeV < pT < 0.8GeV
    histosData.add("fourpion_mass_non_0_charge_domC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // 0.8GeV < pT

    // tpc signal
    histosData.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosData.add("tpcSignal_pions", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

    // tof beta
    histosData.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histosData.add("tofBeta_pions", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

    // Other signals
    histosData.add("FT0A", "T0A amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosData.add("FT0C", "T0C amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosData.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosData.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosData.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});

    // Collin Soper Theta and Phi
    histosData.add("collin_soper_phi_1", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_phi_2", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosData.add("collin_soper_costheta_1", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("collin_soper_costheta_2", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosData.add("phi_vs_costheta_1", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    histosData.add("phi_vs_costheta_2", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});

    // MC Gen Stuff

    // counts
    histosMCgen.add("rhoPrimeCounts", "Total Rho prime Events; Events", kTH1F, {{10, 0, 10}});

    // Track Stuff
    histosMCgen.add("pion_pT", "Generated pT; pT [GeV/c]; Events", kTH1F, {{1000, 0, 10}});
    histosMCgen.add("pion_eta", "Generated Pseudorapidity; #eta; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosMCgen.add("pion_rapidity", "Generated Rapidity; y; Events", kTH1F, {rapidityAxis});

    // Pair Invariant Mass
    histosMCgen.add("twoPion_invMass_pair_1", "Invariant Mass Distribution of 2 pions 1 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCgen.add("twoPion_invMass_pair_2", "Invariant Mass Distribution of 2 pions 2 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCgen.add("twoPion_invMass_pair_3", "Invariant Mass Distribution of 2 pions 3 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCgen.add("twoPion_invMass_pair_4", "Invariant Mass Distribution of 2 pions 4 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});

    // Generated Transverse Momentum, Rapidty and Invariant Mass
    histosMCgen.add("rhoPrime_pT", "Generated pT; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosMCgen.add("rhoPrime_eta", "Generated pT; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosMCgen.add("rhoPrime_rapidity", "Generated pT; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosMCgen.add("rhoPrime_invmass", "Generated pT; pT [GeV/c]; Events", kTH1F, {pTAxis});

    histosMCgen.add("fourPion_pT", "Generated pT; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosMCgen.add("fourPion_eta", "Generated Pseudorapidity; #eta; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosMCgen.add("fourPion_rapidity", "Generated Rapidity; y; Events", kTH1F, {rapidityAxis});
    histosMCgen.add("fourPion_invmass", "Invariant Mass of 4-Pions; m(4-pion); Events", kTH1F, {invMassAxis});

    // Collin Soper Theta and Phi
    histosMCgen.add("fourPion_phi_pair_1", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosMCgen.add("fourPion_phi_pair_2", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosMCgen.add("fourPion_costheta_pair_1", "#theta Distribution;cos(#theta); Events", kTH1F, {cosThetaAxis});
    histosMCgen.add("fourPion_costheta_pair_2", "#theta Distribution;cos(#theta); Events", kTH1F, {cosThetaAxis});
    histosMCgen.add("phi_vs_costheta_1", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    histosMCgen.add("phi_vs_costheta_2", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});

    // MC Reco Stuff

    histosMCreco.add("vertexZ", "Vertex Z; Vertex Z [cm]; Counts", kTH1F, {{1000, -20, 20}});
    histosMCreco.add("occupancy", "Occupancy; Occupancy; Counts", kTH1F, {{1500, 0, 1500}});
    histosMCreco.add("dcaXY", "dcaXY; dcaXY [cm]; Counts", kTH1F, {{10000, -5, 5}});
    histosMCreco.add("dcaZ", "dcaZ; dcaZ [cm]; Counts", kTH1F, {{10000, -10, 10}});
    histosMCreco.add("tpcChi2NCl", "TPC Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{200, 0, 200}});
    histosMCreco.add("itsChi2NCl", "ITS Chi2/NCl; Chi2/NCl; Counts", kTH1F, {{200, 0, 200}});
    histosMCreco.add("tpcNClsFindable", "TPC N Cls Findable; N Cls Findable; Counts", kTH1F, {{200, 0, 200}});

    histosMCreco.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histosMCreco.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histosMCreco.add("EventCounts", "Total Events; Events", kTH1F, {{10, 0, 10}});

    // TPC nSigma
    histosMCreco.add("tpcNSigmaPi_all", "TPC nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tpcNSigmaPi_pions", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TPC nSigma of other particles with selected pion tracks
    histosMCreco.add("tpcNSigmaKa_pions", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tpcNSigmaPr_pions", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tpcNSigmaEl_pions", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tpcNSigmaMu_pions", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TOF nSigma
    histosMCreco.add("tofNSigmaPi_all", "TOF nSigma Pion with track selection; Events", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tofNSigmaPi_pions", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // TOF nSigma of other particles with selected pion tracks
    histosMCreco.add("tofNSigmaKa_pions", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tofNSigmaPr_pions", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tofNSigmaEl_pions", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});
    histosMCreco.add("tofNSigmaMu_pions", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH2F, {{1000, -15, 15}, {1000, 0, 10}});

    // Track Transverse Momentum
    histosMCreco.add("pT_track_all", "pT with track selection; pT [GeV/c]; Counts", kTH1F, {pTAxis});
    histosMCreco.add("pT_track_pions", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});
    histosMCreco.add("pT_track_pions_contributed", "pT with track selection and PID selection of Pi which are contributed to selected event; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Track Rapidity
    histosMCreco.add("rapidity_track_all", "Rapidity with track selection; y; Counts", kTH1F, {rapidityAxis});
    histosMCreco.add("rapidity_track_pions", "Rapidity with track selection and PID selection of Pi; y; Events", kTH1F, {rapidityAxis});
    histosMCreco.add("rapidity_track_pions_contributed", "Rapidity with track selection and PID selection of Pi which are contributed to selected event; y; Events", kTH1F, {rapidityAxis});

    // Zero charge Event Transverse Momentum
    histosMCreco.add("fourpion_pT_0_charge", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Non Zero charge Event Transverse Momentum
    histosMCreco.add("fourpion_pT_non_0_charge", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Events", kTH1F, {pTAxis});

    // Rapidity of 0 charge Events
    histosMCreco.add("fourpion_rap_0_charge_domA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosMCreco.add("fourpion_rap_0_charge_domB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});
    histosMCreco.add("fourpion_rap_0_charge_domC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {{1000, -2.5, 2.5}});

    // Rapidity of non 0 charge Events
    histosMCreco.add("fourpion_rap_non_0_charge_domA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Events", kTH1F, {rapidityAxis});
    histosMCreco.add("fourpion_rap_non_0_charge_domB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c$; y; Events", kTH1F, {rapidityAxis});
    histosMCreco.add("fourpion_rap_non_0_charge_domC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Events", kTH1F, {rapidityAxis});

    // Pair Invariant Mass
    histosMCreco.add("twopion_mass_1", "Invariant Mass Distribution of 2 pions 1 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCreco.add("twopion_mass_2", "Invariant Mass Distribution of 2 pions 2 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCreco.add("twopion_mass_3", "Invariant Mass Distribution of 2 pions 3 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});
    histosMCreco.add("twopion_mass_4", "Invariant Mass Distribution of 2 pions 4 ; m(#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{5000, 0, 5}});

    // Invariant Mass of 0 charge events
    histosMCreco.add("fourpion_mass_0_charge_domA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // pT < 0.15GeV
    histosMCreco.add("fourpion_mass_0_charge_domB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis}); // 0.15GeV < pT < 0.8GeV
    histosMCreco.add("fourpion_mass_0_charge_domC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // 0.8GeV < pT

    // Invariant mass of non 0 charge events
    histosMCreco.add("fourpion_mass_non_0_charge_domA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // pT < 0.15GeV
    histosMCreco.add("fourpion_mass_non_0_charge_domB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis}); // 0.15GeV < pT < 0.8GeV
    histosMCreco.add("fourpion_mass_non_0_charge_domC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {invMassAxis});       // 0.8GeV < pT

    // tpc signal
    histosMCreco.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histosMCreco.add("tpcSignal_pions", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

    // tof beta
    histosMCreco.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histosMCreco.add("tofBeta_pions", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

    // Other signals
    histosMCreco.add("FT0A", "T0A amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosMCreco.add("FT0C", "T0C amplitude", kTH1F, {{2000, 0.0, 500.0}});
    histosMCreco.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosMCreco.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histosMCreco.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});

    // Collin Soper Theta and Phi
    histosMCreco.add("collin_soper_phi_1", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosMCreco.add("collin_soper_phi_2", "#phi Distribution; #phi; Events", kTH1F, {phiAxis});
    histosMCreco.add("collin_soper_costheta_1", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosMCreco.add("collin_soper_costheta_2", "#theta Distribution;cos(#theta); Counts", kTH1F, {cosThetaAxis});
    histosMCreco.add("phi_vs_costheta_1", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});
    histosMCreco.add("phi_vs_costheta_2", "Phi vs cosTheta; #phi; cos(#theta)", kTH2F, {phiAxis, cosThetaAxis});

  } // End of init function
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Calculate the Collins-Soper Frame----------------------------------------------------------------------------------------------------------------------------
  double cosThetaCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    //  TVector3 beta = (-1. / fourpion.E()) * fourpion.Vect();
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
  } // End of cosThetaCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  // Calculate Phi in Collins-Soper Frame------------------------------------------------------------------------------------------------------------------------
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
  } // End of phiCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Filter vertexCut = (nabs(o2::aod::collision::posZ) <= vZCut) && (o2::aod::collision::numContrib == numPVContrib);
  Filter fitcuts = o2::aod::udcollision::totalFV0AmplitudeA < fv0Cut && o2::aod::udcollision::totalFT0AmplitudeA < ft0aCut && o2::aod::udcollision::totalFT0AmplitudeC < ft0cCut;
  Filter zdcCuts = (o2::aod::udzdc::energyCommonZNA < zdcCut) && (o2::aod::udzdc::energyCommonZNC < zdcCut);
  Filter bcSelectionCut = (o2::aod::udcollision::sbp == sbpCut) && (o2::aod::udcollision::itsROFb == itsROFbCut) && (o2::aod::udcollision::vtxITSTPC == vtxITSTPCcut) && (o2::aod::udcollision::tfb == tfbCut);
  Filter onlyPVtracks = o2::aod::udtrack::isPVContributor == true;
  Filter occupCut = nabs(o2::aod::udcollision::occupancyInTime) < occupancyCut;
  using UDtracks = soa::Filtered<soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>>;
  using UDCollisions = soa::Filtered<soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>>; //
  using UDCollision = UDCollisions::iterator;
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Process function--------------------------------------------------------------------------------------------------------------------------------------------------
  void processData(UDCollision const& collision, UDtracks const& tracks)
  {

    int gapSide = collision.gapSide();
    std::vector<float> parameters = {pvCut, dcaZcut, dcaXYcut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, pTcut};
    int truegapSide = sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut);
    histosData.fill(HIST("GapSide"), gapSide);
    histosData.fill(HIST("TrueGapSide"), truegapSide);
    histosData.fill(HIST("EventCounts"), 1);
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
      if (trackselector(t0, parameters)) {
        selectedTracks.push_back(t0);
        if (selectionPIDPion(t0, true, nSigmaTPCcut, nSigmaTOFcut)) {
          selectedPionTracks.push_back(t0);
          if (t0.sign() == 1) {
            selectedPionPlusTracks.push_back(t0);
          }
          if (t0.sign() == -1) {
            selectedPionMinusTracks.push_back(t0);
          }
        } // End of Selection PID Pion
      } // End of track selections
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
      histosData.fill(HIST("rapidity_track_all"), selectedTrackVector.Rapidity());
      histosData.fill(HIST("itsChi2NCl"), selectedTracks[i].itsChi2NCl());
      histosData.fill(HIST("tpcChi2NCl"), selectedTracks[i].tpcChi2NCl());
      histosData.fill(HIST("tpcNClsFindable"), selectedTracks[i].tpcNClsFindable());
      histosData.fill(HIST("dcaXY"), selectedTracks[i].dcaXY());
      histosData.fill(HIST("dcaZ"), selectedTracks[i].dcaZ());
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
      histosData.fill(HIST("rapidity_track_pions"), selectedPionTrackVector.Rapidity());
    } // End of loop over tracks with selection and PID selection of Pions

    if (numSelectedPionTracks != numFourPionTracks) {
      return;
    }

    histosData.fill(HIST("EventCounts"), 2);

    // Check if there is at least one track with TOF in the selected events, otherwise return
    bool hasAtleastOneTOF = false;
    for (int i = 0; i < numPiPlusTracks; i++) {
      if (selectedPionPlusTracks[i].hasTOF() == true) {
        hasAtleastOneTOF = true;
        break;
      }
    }
    if ((!hasAtleastOneTOF) && (checkOneTof == 1)) {
      return;
    }
    histosData.fill(HIST("EventCounts"), 3);

    // Selecting Events with net charge = 0
    if (numPionMinusTracks == numPiMinus && numPiPlusTracks == numPiPlus) {

      histosData.fill(HIST("EventCounts"), 4);

      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      ROOT::Math::PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

      histosData.fill(HIST("pT_track_pions_contributed"), p1.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p2.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p3.Pt());
      histosData.fill(HIST("pT_track_pions_contributed"), p4.Pt());

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

      double fourPiPhiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
      double fourPiPhiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
      double fourPiCosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
      double fourPiCosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

      sigFromData(
        collision.posX(), collision.posY(), collision.posZ(),
        collision.totalFV0AmplitudeA(), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
        collision.timeFV0A(), collision.timeFT0A(), collision.timeFT0C(), collision.timeFDDA(), collision.timeFDDC(), collision.timeZNA(), collision.timeZNC(), collision.occupancyInTime(),
        selectedPionPlusTracks[0].dcaXY(), selectedPionPlusTracks[1].dcaXY(), selectedPionMinusTracks[0].dcaXY(), selectedPionMinusTracks[1].dcaXY(),

        selectedPionPlusTracks[0].dcaZ(), selectedPionPlusTracks[1].dcaZ(), selectedPionMinusTracks[0].dcaZ(), selectedPionMinusTracks[1].dcaZ(),

        selectedPionPlusTracks[0].tpcNSigmaPi(), selectedPionPlusTracks[1].tpcNSigmaPi(), selectedPionMinusTracks[0].tpcNSigmaPi(), selectedPionMinusTracks[1].tpcNSigmaPi(),

        selectedPionPlusTracks[0].tpcNSigmaKa(), selectedPionPlusTracks[1].tpcNSigmaKa(), selectedPionMinusTracks[0].tpcNSigmaKa(), selectedPionMinusTracks[1].tpcNSigmaKa(),

        selectedPionPlusTracks[0].tpcNSigmaPr(), selectedPionPlusTracks[1].tpcNSigmaPr(), selectedPionMinusTracks[0].tpcNSigmaPr(), selectedPionMinusTracks[1].tpcNSigmaPr(),

        selectedPionPlusTracks[0].tpcNSigmaEl(), selectedPionPlusTracks[1].tpcNSigmaEl(), selectedPionMinusTracks[0].tpcNSigmaEl(), selectedPionMinusTracks[1].tpcNSigmaEl(),

        selectedPionPlusTracks[0].tpcNSigmaMu(), selectedPionPlusTracks[1].tpcNSigmaMu(), selectedPionMinusTracks[0].tpcNSigmaMu(), selectedPionMinusTracks[1].tpcNSigmaMu(),

        selectedPionPlusTracks[0].tpcChi2NCl(), selectedPionPlusTracks[1].tpcChi2NCl(), selectedPionMinusTracks[0].tpcChi2NCl(), selectedPionMinusTracks[1].tpcChi2NCl(),

        selectedPionPlusTracks[0].tpcNClsFindable(), selectedPionPlusTracks[1].tpcNClsFindable(), selectedPionMinusTracks[0].tpcNClsFindable(), selectedPionMinusTracks[1].tpcNClsFindable(),

        selectedPionPlusTracks[0].itsChi2NCl(), selectedPionPlusTracks[1].itsChi2NCl(), selectedPionMinusTracks[0].itsChi2NCl(), selectedPionMinusTracks[1].itsChi2NCl(),

        p1.Pt(), p2.Pt(), p3.Pt(), p4.Pt(),
        p1.Eta(), p2.Eta(), p3.Eta(), p4.Eta(),
        p1.Phi(), p2.Phi(), p3.Phi(), p4.Phi(),
        p1.Rapidity(), p2.Rapidity(), p3.Rapidity(), p4.Rapidity(),
        p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(),
        fourPiPhiPair1, fourPiPhiPair2, fourPiCosThetaPair1, fourPiCosThetaPair2);

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_0_charge"), p1234.Pt());
        if (p1234.Pt() < rhoPtCut) {

          histosData.fill(HIST("EventCounts"), 5);

          histosData.fill(HIST("fourpion_rap_0_charge_domA"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_0_charge_domA"), p1234.M());

          histosData.fill(HIST("twopion_mass_1"), (p1 + p3).M());
          histosData.fill(HIST("twopion_mass_2"), (p1 + p4).M());
          histosData.fill(HIST("twopion_mass_3"), (p2 + p3).M());
          histosData.fill(HIST("twopion_mass_4"), (p2 + p4).M());

          histosData.fill(HIST("collin_soper_phi_1"), fourPiPhiPair1);
          histosData.fill(HIST("collin_soper_phi_2"), fourPiPhiPair2);
          histosData.fill(HIST("collin_soper_costheta_1"), fourPiCosThetaPair1);
          histosData.fill(HIST("collin_soper_costheta_2"), fourPiCosThetaPair2);

          histosData.fill(HIST("phi_vs_costheta_1"), fourPiPhiPair1, fourPiCosThetaPair1);
          histosData.fill(HIST("phi_vs_costheta_2"), fourPiPhiPair2, fourPiCosThetaPair2);
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosData.fill(HIST("fourpion_rap_0_charge_domB"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosData.fill(HIST("fourpion_rap_0_charge_domC"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTracks != numPiMinus && numPiPlusTracks != numPiPlus) {

      histosData.fill(HIST("EventCounts"), 6);

      ROOT::Math::PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);

      ROOT::Math::PxPyPzMVector p1234 = p1 + p2 + p3 + p4;

      bkgFromData(
        collision.posX(), collision.posY(), collision.posZ(),
        collision.totalFV0AmplitudeA(), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
        collision.timeFV0A(), collision.timeFT0A(), collision.timeFT0C(), collision.timeFDDA(), collision.timeFDDC(),
        collision.timeZNA(), collision.timeZNC(), collision.occupancyInTime(),
        selectedPionTracks[0].dcaXY(), selectedPionTracks[1].dcaXY(), selectedPionTracks[0].dcaXY(), selectedPionTracks[1].dcaXY(),
        selectedPionTracks[0].dcaZ(), selectedPionTracks[1].dcaZ(), selectedPionTracks[0].dcaZ(), selectedPionTracks[1].dcaZ(),
        selectedPionTracks[0].tpcNSigmaPi(), selectedPionTracks[1].tpcNSigmaPi(), selectedPionTracks[0].tpcNSigmaPi(), selectedPionTracks[1].tpcNSigmaPi(),
        selectedPionTracks[0].tpcNSigmaKa(), selectedPionTracks[1].tpcNSigmaKa(), selectedPionTracks[0].tpcNSigmaKa(), selectedPionTracks[1].tpcNSigmaKa(),
        selectedPionTracks[0].tpcNSigmaPr(), selectedPionTracks[1].tpcNSigmaPr(), selectedPionTracks[0].tpcNSigmaPr(), selectedPionTracks[1].tpcNSigmaPr(),
        selectedPionTracks[0].tpcNSigmaEl(), selectedPionTracks[1].tpcNSigmaEl(), selectedPionTracks[0].tpcNSigmaEl(), selectedPionTracks[1].tpcNSigmaEl(),
        selectedPionTracks[0].tpcNSigmaMu(), selectedPionTracks[1].tpcNSigmaMu(), selectedPionTracks[0].tpcNSigmaMu(), selectedPionTracks[1].tpcNSigmaMu(),
        selectedPionTracks[0].itsChi2NCl(), selectedPionTracks[1].itsChi2NCl(), selectedPionTracks[0].itsChi2NCl(), selectedPionTracks[1].itsChi2NCl(),
        selectedPionTracks[0].tpcChi2NCl(), selectedPionTracks[1].tpcChi2NCl(), selectedPionTracks[0].tpcChi2NCl(), selectedPionTracks[1].tpcChi2NCl(),
        selectedPionTracks[0].tpcNClsFindable(), selectedPionTracks[1].tpcNClsFindable(), selectedPionTracks[0].tpcNClsFindable(), selectedPionTracks[1].tpcNClsFindable(),
        p1.Pt(), p2.Pt(), p3.Pt(), p4.Pt(),
        p1.Eta(), p2.Eta(), p3.Eta(), p4.Eta(),
        p1.Phi(), p2.Phi(), p3.Phi(), p4.Phi(),
        p1.Rapidity(), p2.Rapidity(), p3.Rapidity(), p4.Rapidity(),
        p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M());

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosData.fill(HIST("fourpion_pT_non_0_charge"), p1234.Pt());

        if (p1234.Pt() < rhoPtCut) {
          histosData.fill(HIST("EventCounts"), 7);
          histosData.fill(HIST("fourpion_rap_non_0_charge_domA"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_non_0_charge_domA"), p1234.M());
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosData.fill(HIST("fourpion_rap_non_0_charge_domB"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_non_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosData.fill(HIST("fourpion_rap_non_0_charge_domC"), p1234.Rapidity());
          histosData.fill(HIST("fourpion_mass_non_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection
    } // End of Analysis for non 0 charge events
  } // End of 4 Pion Analysis Process function for Data
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processData, "The Process for 4 Pion Analysis from data", false);
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of MC Generation function-----------------------------------------------------------------------------------------------------------------------------------------------
  void processMCgen(aod::UDMcCollisions::iterator const&, aod::UDMcParticles const& mcParts)
  {
    std::vector<ROOT::Math::PxPyPzMVector> piPlusvectors;
    std::vector<ROOT::Math::PxPyPzMVector> piMinusvectors;
    TVector3 particleVector;

    for (const auto& particle : mcParts) {

      if ((particle.pdgCode() != rhoPrime) || (particle.daughters_as<aod::UDMcParticles>().size() != numFourPionTracks)) {
        continue;
      }

      particleVector.SetXYZ(particle.px(), particle.py(), particle.pz());
      histosMCgen.fill(HIST("rhoPrimeCounts"), 1);
      histosMCgen.fill(HIST("rhoPrime_pT"), particleVector.Pt());
      histosMCgen.fill(HIST("rhoPrime_eta"), particleVector.Eta());

      for (const auto& daughter : particle.daughters_as<aod::UDMcParticles>()) {
        ROOT::Math::PxPyPzMVector daughterVector(daughter.px(), daughter.py(), daughter.pz(), o2::constants::physics::MassPionCharged);
        if (daughter.pdgCode() == PDG_t::kPiPlus) {
          piPlusvectors.push_back(daughterVector);
        }
        if (daughter.pdgCode() == PDG_t::kPiMinus) {
          piMinusvectors.push_back(daughterVector);
        }
      } // End of loop over daughters

    } // End of loop over MC particles

    if (static_cast<int>(piPlusvectors.size()) != numPiPlus || static_cast<int>(piMinusvectors.size()) != numPiMinus) {
      return;
    }

    ROOT::Math::PxPyPzMVector p1234 = piPlusvectors[0] + piPlusvectors[1] + piMinusvectors[0] + piMinusvectors[1];

    histosMCgen.fill(HIST("pion_pT"), piPlusvectors[0].Pt());
    histosMCgen.fill(HIST("pion_pT"), piPlusvectors[1].Pt());
    histosMCgen.fill(HIST("pion_pT"), piMinusvectors[0].Pt());
    histosMCgen.fill(HIST("pion_pT"), piMinusvectors[1].Pt());

    histosMCgen.fill(HIST("pion_eta"), piPlusvectors[0].Eta());
    histosMCgen.fill(HIST("pion_eta"), piPlusvectors[1].Eta());
    histosMCgen.fill(HIST("pion_eta"), piMinusvectors[0].Eta());
    histosMCgen.fill(HIST("pion_eta"), piMinusvectors[1].Eta());

    histosMCgen.fill(HIST("pion_rapidity"), piPlusvectors[0].Rapidity());
    histosMCgen.fill(HIST("pion_rapidity"), piPlusvectors[1].Rapidity());
    histosMCgen.fill(HIST("pion_rapidity"), piMinusvectors[0].Rapidity());
    histosMCgen.fill(HIST("pion_rapidity"), piMinusvectors[1].Rapidity());

    histosMCgen.fill(HIST("fourPion_pT"), p1234.Pt());
    histosMCgen.fill(HIST("fourPion_eta"), p1234.Eta());
    histosMCgen.fill(HIST("fourPion_rapidity"), p1234.Rapidity());
    histosMCgen.fill(HIST("fourPion_invmass"), p1234.M());

    histosMCgen.fill(HIST("twoPion_invMass_pair_1"), (piPlusvectors[0] + piMinusvectors[0]).M());
    histosMCgen.fill(HIST("twoPion_invMass_pair_2"), (piPlusvectors[0] + piMinusvectors[1]).M());
    histosMCgen.fill(HIST("twoPion_invMass_pair_3"), (piPlusvectors[1] + piMinusvectors[0]).M());
    histosMCgen.fill(HIST("twoPion_invMass_pair_4"), (piPlusvectors[1] + piMinusvectors[1]).M());

    ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

    k1.SetCoordinates(piPlusvectors[0].Pt(), piPlusvectors[0].Eta(), piPlusvectors[0].Phi(), o2::constants::physics::MassPionCharged);
    k2.SetCoordinates(piPlusvectors[1].Pt(), piPlusvectors[1].Eta(), piPlusvectors[1].Phi(), o2::constants::physics::MassPionCharged);
    k3.SetCoordinates(piMinusvectors[0].Pt(), piMinusvectors[0].Eta(), piMinusvectors[0].Phi(), o2::constants::physics::MassPionCharged);
    k4.SetCoordinates(piMinusvectors[1].Pt(), piMinusvectors[1].Eta(), piMinusvectors[1].Phi(), o2::constants::physics::MassPionCharged);

    k1234 = k1 + k2 + k3 + k4;

    k13 = k1 + k3;
    k14 = k1 + k4;
    k23 = k2 + k3;
    k24 = k2 + k4;

    auto phiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
    auto phiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
    auto cosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
    auto cosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

    generatedMC(
      piPlusvectors[0].Pt(), piPlusvectors[1].Pt(), piMinusvectors[0].Pt(), piMinusvectors[1].Pt(),
      piPlusvectors[0].Eta(), piPlusvectors[1].Eta(), piMinusvectors[0].Eta(), piMinusvectors[1].Eta(),
      piPlusvectors[0].Phi(), piPlusvectors[1].Phi(), piMinusvectors[0].Phi(), piMinusvectors[1].Phi(),
      piPlusvectors[0].Rapidity(), piPlusvectors[1].Rapidity(), piMinusvectors[0].Rapidity(), piMinusvectors[1].Rapidity(),
      p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(),
      phiPair1, phiPair2, cosThetaPair1, cosThetaPair2);

    histosMCgen.fill(HIST("fourPion_phi_pair_1"), phiPair1);
    histosMCgen.fill(HIST("fourPion_phi_pair_2"), phiPair2);
    histosMCgen.fill(HIST("fourPion_costheta_pair_1"), cosThetaPair1);
    histosMCgen.fill(HIST("fourPion_costheta_pair_2"), cosThetaPair2);
    histosMCgen.fill(HIST("phi_vs_costheta_1"), phiPair1, cosThetaPair1);
    histosMCgen.fill(HIST("phi_vs_costheta_2"), phiPair2, cosThetaPair2);

  } // End of 4 Pion MC Generation Process function
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processMCgen, "The Process for 4 Pion Analysis from MC Generation", true);

  // Begin of MC Reconstruction function-----------------------------------------------------------------------------------------------------------------------------------------------
  using CollisionStuff = soa::Join<aod::UDCollisions_001, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDMcCollsLabels>;
  using CollisionTotal = CollisionStuff::iterator;
  using TrackStuff = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;

  void processMCrec(CollisionTotal const& collision, TrackStuff const& tracks)
  {

    if (!collision.has_udMcCollision()) {
      return;
    }

    int gapSide = collision.gapSide();
    std::vector<float> parameters = {pvCut, dcaZcut, dcaXYcut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, pTcut};
    int truegapSide = sgSelector.trueGap(collision, fv0Cut, ft0aCut, ft0cCut, zdcCut);
    histosMCreco.fill(HIST("GapSide"), gapSide);
    histosMCreco.fill(HIST("TrueGapSide"), truegapSide);
    histosMCreco.fill(HIST("EventCounts"), 1);
    histosMCreco.fill(HIST("vertexZ"), collision.posZ());
    // histosMCreco.fill(HIST("occupancy"), collision.occupancyInTime());
    histosMCreco.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histosMCreco.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histosMCreco.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histosMCreco.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histosMCreco.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    std::vector<decltype(tracks.begin())> selectedTracks;
    std::vector<decltype(tracks.begin())> selectedPionTracks;
    std::vector<decltype(tracks.begin())> selectedPionPlusTracks;
    std::vector<decltype(tracks.begin())> selectedPionMinusTracks;

    for (const auto& t0 : tracks) {
      if (trackselector(t0, parameters) && t0.has_udMcParticle()) {
        selectedTracks.push_back(t0);
        if (selectionPIDPion(t0, true, nSigmaTPCcut, nSigmaTOFcut)) {
          selectedPionTracks.push_back(t0);
          if (t0.sign() == 1) {
            selectedPionPlusTracks.push_back(t0);
          }
          if (t0.sign() == -1) {
            selectedPionMinusTracks.push_back(t0);
          }
        } // End of Selection PID Pion
      } // End of track selections
    } // End of loop over tracks

    int numSelectedTracks = static_cast<int>(selectedTracks.size());
    int numSelectedPionTracks = static_cast<int>(selectedPionTracks.size());
    int numPiPlusTracks = static_cast<int>(selectedPionPlusTracks.size());
    int numPionMinusTRacks = static_cast<int>(selectedPionMinusTracks.size());

    for (int i = 0; i < numSelectedTracks; i++) {
      ROOT::Math::PxPyPzMVector selectedTrackVector(selectedTracks[i].px(), selectedTracks[i].py(), selectedTracks[i].pz(), o2::constants::physics::MassPionCharged);
      histosMCreco.fill(HIST("tpcSignal"), selectedTrackVector.P(), selectedTracks[i].tpcSignal());
      histosMCreco.fill(HIST("tofBeta"), selectedTrackVector.P(), selectedTracks[i].beta());
      histosMCreco.fill(HIST("tpcNSigmaPi_all"), selectedTracks[i].tpcNSigmaPi(), selectedTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaPi_all"), selectedTracks[i].tofNSigmaPi(), selectedTrackVector.Pt());
      histosMCreco.fill(HIST("pT_track_all"), selectedTrackVector.Pt());
      histosMCreco.fill(HIST("rapidity_track_all"), selectedTrackVector.Rapidity());
      histosMCreco.fill(HIST("itsChi2NCl"), selectedTracks[i].itsChi2NCl());
      histosMCreco.fill(HIST("tpcChi2NCl"), selectedTracks[i].tpcChi2NCl());
      histosMCreco.fill(HIST("tpcNClsFindable"), selectedTracks[i].tpcNClsFindable());
      histosMCreco.fill(HIST("dcaXY"), selectedTracks[i].dcaXY());
      histosMCreco.fill(HIST("dcaZ"), selectedTracks[i].dcaZ());
    } // End of loop over tracks with selection only

    for (int i = 0; i < numSelectedPionTracks; i++) {
      ROOT::Math::PxPyPzMVector selectedPionTrackVector(selectedPionTracks[i].px(), selectedPionTracks[i].py(), selectedPionTracks[i].pz(), o2::constants::physics::MassPionCharged);
      histosMCreco.fill(HIST("tpcSignal_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].tpcSignal());
      histosMCreco.fill(HIST("tofBeta_pions"), selectedPionTrackVector.P(), selectedPionTracks[i].beta());
      histosMCreco.fill(HIST("tpcNSigmaPi_pions"), selectedPionTracks[i].tpcNSigmaPi(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tpcNSigmaKa_pions"), selectedPionTracks[i].tpcNSigmaKa(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tpcNSigmaPr_pions"), selectedPionTracks[i].tpcNSigmaPr(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tpcNSigmaEl_pions"), selectedPionTracks[i].tpcNSigmaEl(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tpcNSigmaMu_pions"), selectedPionTracks[i].tpcNSigmaMu(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaPi_pions"), selectedPionTracks[i].tofNSigmaPi(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaKa_pions"), selectedPionTracks[i].tofNSigmaKa(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaPr_pions"), selectedPionTracks[i].tofNSigmaPr(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaEl_pions"), selectedPionTracks[i].tofNSigmaEl(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("tofNSigmaMu_pions"), selectedPionTracks[i].tofNSigmaMu(), selectedPionTrackVector.Pt());
      histosMCreco.fill(HIST("pT_track_pions"), std::sqrt(selectedPionTracks[i].px() * selectedPionTracks[i].px() + selectedPionTracks[i].py() * selectedPionTracks[i].py()));
      histosMCreco.fill(HIST("rapidity_track_pions"), selectedPionTrackVector.Rapidity());
    } // End of loop over tracks with selection and PID selection of Pions

    if (numSelectedPionTracks != numFourPionTracks) {
      return;
    }

    // Check if there is at least one track with TOF in the selected events, otherwise return
    bool hasAtleastOneTOF = false;
    for (int i = 0; i < numPiPlusTracks; i++) {
      if (selectedPionPlusTracks[i].hasTOF() == true) {
        hasAtleastOneTOF = true;
        break;
      }
    }
    if ((!hasAtleastOneTOF) && (checkOneTof == 1)) {
      return;
    }
    histosMCreco.fill(HIST("EventCounts"), 3);

    // Selecting Events with net charge = 0
    if (numPionMinusTRacks == numPiMinus && numPiPlusTracks == numPiPlus) {

      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      ROOT::Math::PxPyPzMVector p1(selectedPionPlusTracks[0].px(), selectedPionPlusTracks[0].py(), selectedPionPlusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionPlusTracks[1].px(), selectedPionPlusTracks[1].py(), selectedPionPlusTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionMinusTracks[0].px(), selectedPionMinusTracks[0].py(), selectedPionMinusTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionMinusTracks[1].px(), selectedPionMinusTracks[1].py(), selectedPionMinusTracks[1].pz(), o2::constants::physics::MassPionCharged);

      histosMCreco.fill(HIST("pT_track_pions_contributed"), p1.Pt());
      histosMCreco.fill(HIST("pT_track_pions_contributed"), p2.Pt());
      histosMCreco.fill(HIST("pT_track_pions_contributed"), p3.Pt());
      histosMCreco.fill(HIST("pT_track_pions_contributed"), p4.Pt());

      histosMCreco.fill(HIST("rapidity_track_pions_contributed"), p1.Rapidity());
      histosMCreco.fill(HIST("rapidity_track_pions_contributed"), p2.Rapidity());
      histosMCreco.fill(HIST("rapidity_track_pions_contributed"), p3.Rapidity());
      histosMCreco.fill(HIST("rapidity_track_pions_contributed"), p4.Rapidity());

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

      double phiPair1 = phiCollinsSoperFrame(k13, k24, k1234);
      double phiPair2 = phiCollinsSoperFrame(k14, k23, k1234);
      double cosThetaPair1 = cosThetaCollinsSoperFrame(k13, k24, k1234);
      double cosThetaPair2 = cosThetaCollinsSoperFrame(k14, k23, k1234);

      sigFromMC(
        collision.posX(), collision.posY(), collision.posZ(),
        collision.totalFV0AmplitudeA(), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(),
        collision.timeFV0A(), collision.timeFT0A(), collision.timeFT0C(), collision.timeFDDA(), collision.timeFDDC(),
        collision.timeZNA(), collision.timeZNC(), 100,
        selectedPionPlusTracks[0].dcaXY(), selectedPionPlusTracks[1].dcaXY(), selectedPionMinusTracks[0].dcaXY(), selectedPionMinusTracks[1].dcaXY(),
        selectedPionPlusTracks[0].dcaZ(), selectedPionPlusTracks[1].dcaZ(), selectedPionMinusTracks[0].dcaZ(), selectedPionMinusTracks[1].dcaZ(),
        selectedPionPlusTracks[0].tpcNSigmaPi(), selectedPionPlusTracks[1].tpcNSigmaPi(), selectedPionMinusTracks[0].tpcNSigmaPi(), selectedPionMinusTracks[1].tpcNSigmaPi(),
        selectedPionPlusTracks[0].tpcNSigmaKa(), selectedPionPlusTracks[1].tpcNSigmaKa(), selectedPionMinusTracks[0].tpcNSigmaKa(), selectedPionMinusTracks[1].tpcNSigmaKa(),
        selectedPionPlusTracks[0].tpcNSigmaPr(), selectedPionPlusTracks[1].tpcNSigmaPr(), selectedPionMinusTracks[0].tpcNSigmaPr(), selectedPionMinusTracks[1].tpcNSigmaPr(),
        selectedPionPlusTracks[0].tpcNSigmaEl(), selectedPionPlusTracks[1].tpcNSigmaEl(), selectedPionMinusTracks[0].tpcNSigmaEl(), selectedPionMinusTracks[1].tpcNSigmaEl(),
        selectedPionPlusTracks[0].tpcNSigmaMu(), selectedPionPlusTracks[1].tpcNSigmaMu(), selectedPionMinusTracks[0].tpcNSigmaMu(), selectedPionMinusTracks[1].tpcNSigmaMu(),
        selectedPionPlusTracks[0].tpcChi2NCl(), selectedPionPlusTracks[1].tpcChi2NCl(), selectedPionMinusTracks[0].tpcChi2NCl(), selectedPionMinusTracks[1].tpcChi2NCl(),
        selectedPionPlusTracks[0].tpcNClsFindable(), selectedPionPlusTracks[1].tpcNClsFindable(), selectedPionMinusTracks[0].tpcNClsFindable(), selectedPionMinusTracks[1].tpcNClsFindable(),
        selectedPionPlusTracks[0].itsChi2NCl(), selectedPionPlusTracks[1].itsChi2NCl(), selectedPionMinusTracks[0].itsChi2NCl(), selectedPionMinusTracks[1].itsChi2NCl(),
        p1.Pt(), p2.Pt(), p3.Pt(), p4.Pt(),
        p1.Eta(), p2.Eta(), p3.Eta(), p4.Eta(),
        p1.Phi(), p2.Phi(), p3.Phi(), p4.Phi(),
        p1.Rapidity(), p2.Rapidity(), p3.Rapidity(), p4.Rapidity(),
        p1234.Pt(), p1234.Eta(), p1234.Phi(), p1234.Rapidity(), p1234.M(),
        phiPair1, phiPair2, cosThetaPair1, cosThetaPair2);

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosMCreco.fill(HIST("fourpion_pT_0_charge"), p1234.Pt());
        if (p1234.Pt() < rhoPtCut) {
          histosMCreco.fill(HIST("fourpion_rap_0_charge_domA"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_0_charge_domA"), p1234.M());

          histosMCreco.fill(HIST("twopion_mass_1"), (p1 + p3).M());
          histosMCreco.fill(HIST("twopion_mass_2"), (p1 + p4).M());
          histosMCreco.fill(HIST("twopion_mass_3"), (p2 + p3).M());
          histosMCreco.fill(HIST("twopion_mass_4"), (p2 + p4).M());

          histosMCreco.fill(HIST("collin_soper_phi_1"), phiPair1);
          histosMCreco.fill(HIST("collin_soper_phi_2"), phiPair2);
          histosMCreco.fill(HIST("collin_soper_costheta_1"), cosThetaPair1);
          histosMCreco.fill(HIST("collin_soper_costheta_2"), cosThetaPair2);
          histosMCreco.fill(HIST("phi_vs_costheta_1"), phiPair1, cosThetaPair1);
          histosMCreco.fill(HIST("phi_vs_costheta_2"), phiPair2, cosThetaPair2);
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosMCreco.fill(HIST("fourpion_rap_0_charge_domB"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosMCreco.fill(HIST("fourpion_rap_0_charge_domC"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (numPionMinusTRacks != numPiMinus && numPiPlusTracks != numPiPlus) {
      ROOT::Math::PxPyPzMVector p1(selectedPionTracks[0].px(), selectedPionTracks[0].py(), selectedPionTracks[0].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p2(selectedPionTracks[1].px(), selectedPionTracks[1].py(), selectedPionTracks[1].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p3(selectedPionTracks[2].px(), selectedPionTracks[2].py(), selectedPionTracks[2].pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector p4(selectedPionTracks[3].px(), selectedPionTracks[3].py(), selectedPionTracks[3].pz(), o2::constants::physics::MassPionCharged);

      ROOT::Math::PxPyPzMVector p1234 = p1 + p2 + p3 + p4;

      if (std::fabs(p1234.Rapidity()) < rhoRapCut) {
        histosMCreco.fill(HIST("fourpion_pT_non_0_charge"), p1234.Pt());

        if (p1234.Pt() < rhoPtCut) {
          histosMCreco.fill(HIST("fourpion_rap_non_0_charge_domA"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_non_0_charge_domA"), p1234.M());
        }
        if (p1234.Pt() > rhoPtCut && p1234.Pt() < zeroPointEight) {
          histosMCreco.fill(HIST("fourpion_rap_non_0_charge_domB"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_non_0_charge_domB"), p1234.M());
        }
        if (p1234.Pt() > zeroPointEight) {
          histosMCreco.fill(HIST("fourpion_rap_non_0_charge_domC"), p1234.Rapidity());
          histosMCreco.fill(HIST("fourpion_mass_non_0_charge_domC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for non 0 charge events

  } // End of 4 Pion Analysis Process function for MC Reconstruction
  PROCESS_SWITCH(ExclusiveRhoTo4Pi, processMCrec, "The Process for 4 Pion Analysis from MC Reconstruction", false);

}; // End of Struct exclusiveRhoTo4Pi
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ExclusiveRhoTo4Pi>(cfgc)};
}
