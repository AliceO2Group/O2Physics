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
// \Single Gap Event Analyzer
// \author Anantha Padmanabhan M Nair, anantha.manoj.nair@cern.ch
// \since  May 2024

#include <cstdlib>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "Common/DataModel/PIDResponse.h"
#include <TString.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct exclusiveRhoTo4Pi {
  SGSelector sgSelector;
  HistogramRegistry histos{"HistoReg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 150., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};

  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 3.2, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 2.4, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 80, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0, "Track Pt"};

  Configurable<float> nSigmaTPC_cut{"nsigmatpccut", 3, "TPC cut"};
  Configurable<float> nSigmaTOF_cut{"nsigmatofcut", 3, "TOF cut"};
  Configurable<bool> StrictEventSelection{"StrictEventSelection", true, "Event Selection"};

  Configurable<int> nBins_pT{"nBins_pT", 1000, "Number of bins for pT"};
  Configurable<int> nBins_IM{"nBins_IM", 1000, "Number of bins for Invariant Mass"};
  Configurable<int> nBins_y{"nBins_y", 1000, "Number of bins for Rapidity"};
  Configurable<int> nBins_Phi{"nBins_phi", 360, "Number of bins for Phi"};
  Configurable<int> nBins_CosTheta{"nBins_cosTheta", 360, "Number of bins for cos Theta"};
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Init Function-----------------------------------------------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {

    histos.add("GapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histos.add("TrueGapSide", "Gap Side; Events", kTH1F, {{4, -1.5, 2.5}});
    histos.add("EventCounts", "Total Events; Events", kTH1F, {{10, 0, 10}}); // 2=#Events, 3= #Selected Events, 5=#Selected Events with Zero Net charge, 6=Selected Events with Non-Zero Net charge

    // TPC nSigma
    histos.add("tpcNSigmaPi_WOTS", "TPC nSigma Pion without track selection; Events", kTH1F, {{1000, -15, 15}});
    histos.add("tpcNSigmaPi_WTS", "TPC nSigma Pion with track selection; Events", kTH1F, {{1000, -15, 15}});
    histos.add("tpcNSigmaPi_WTS_PID_Pi", "TPC nSigma Pion with track selection and PID Selection of Pi; Entries", kTH1F, {{1000, -15, 15}});

    // TPC nSigma of other particles with selected pion tracks
    histos.add("tpcNSigmaPi_WTS_PID_Pi_Ka", "TPC nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tpcNSigmaPi_WTS_PID_Pi_Pr", "TPC nSigma Proton with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tpcNSigmaPi_WTS_PID_Pi_El", "TPC nSigma Electron with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tpcNSigmaPi_WTS_PID_Pi_Mu", "TPC nSigma Muon with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});

    // TOF nSigma
    histos.add("tofNSigmaPi_WTS", "TOF nSigma Pion with track selection; Events", kTH1F, {{1000, -15, 15}});
    histos.add("tofNSigmaPi_WOTS", "TOF nSigma Pion without track selection; Events", kTH1F, {{1000, -15, 15}});
    histos.add("tofNSigmaPi_WTS_PID_Pi", "TOF nSigma Pion with track selection and PID Selection of Pi; Entries", kTH1F, {{1000, -15, 15}});

    // TOF nSigma of other particles with selected pion tracks
    histos.add("tofNSigmaPi_WTS_PID_Pi_Ka", "TOF nSigma Kaon with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tofNSigmaPi_WTS_PID_Pi_Pr", "TOF nSigma Proton with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tofNSigmaPi_WTS_PID_Pi_El", "TOF nSigma Electron with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});
    histos.add("tofNSigmaPi_WTS_PID_Pi_Mu", "TOF nSigma Muon with track selection and PID Selection of Pion; Entries", kTH1F, {{1000, -15, 15}});

    // Track Transverse Momentum
    histos.add("pT_track_WOTS", "pT without track selection; pT [GeV/c]; Events", kTH1F, {{nBins_pT, 0, 2}});
    histos.add("pT_track_WTS", "pT with track selection; pT [GeV/c]; Events", kTH1F, {{nBins_pT, 0, 2}});
    histos.add("pT_track_WTS_PID_Pi", "pT with track selection and PID selection of Pi; pT [GeV/c]; Events", kTH1F, {{nBins_pT, 0, 2}});

    // Zero charge Event Transverse Momentum
    histos.add("pT_event_0charge_WTS_PID_Pi", "Event pT in 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Counts", kTH1F, {{nBins_pT, 0, 2}});

    // Non Zero charge Event Transverse Momentum
    histos.add("pT_event_non0charge_WTS_PID_Pi", "Event pT in Non 0 Charge Events With Track Selection and PID Selection of Pi; pT [GeV/c]; Counts", kTH1F, {{nBins_pT, 0, 2}});

    // Rapidity of 0 charge Events
    histos.add("rapidity_event_0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Counts", kTH1F, {{1000, -2.5, 2.5}});
    histos.add("rapidity_event_0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; y; Counts", kTH1F, {{1000, -2.5, 2.5}});
    histos.add("rapidity_event_0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Counts", kTH1F, {{1000, -2.5, 2.5}});

    // Rapidity of non 0 charge Events
    histos.add("rapidity_event_non0charge_WTS_PID_Pi_domainA", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} < 0.15 GeV/c; y; Counts", kTH1F, {{nBins_y, -2.5, 2.5}});
    histos.add("rapidity_event_non0charge_WTS_PID_Pi_domainB", "Rapidity of Events With Track Selection and PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c$; y; Counts", kTH1F, {{nBins_y, -2.5, 2.5}});
    histos.add("rapidity_event_non0charge_WTS_PID_Pi_domainC", "Rapidity of Events With Track Selection and PID Selection of Pi for p_{T} > 0.80 GeV/c; y; Counts", kTH1F, {{nBins_y, -2.5, 2.5}});

    // Invariant Mass of 0 charge events
    histos.add("invMass_event_0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}});       // pT < 0.15GeV
    histos.add("invMass_event_0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
    histos.add("invMass_event_0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}});       // 0.8GeV < pT

    // Invariant mass of non 0 charge events
    histos.add("invMass_event_non0charge_WTS_PID_Pi_domainA", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} < 0.15 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}});       // pT < 0.15GeV
    histos.add("invMass_event_non0charge_WTS_PID_Pi_domainB", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for 0.15< p_{T} < 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}}); // 0.15GeV < pT < 0.8GeV
    histos.add("invMass_event_non0charge_WTS_PID_Pi_domainC", "Invariant Mass Distribution of non 0 charge Events with PID Selection of Pi for p_{T} > 0.80 GeV/c; m(#pi^{+}#pi^{-}#pi^{+}#pi^{-}) [GeV/c]", kTH1F, {{nBins_IM, 0.8, 2.5}});       // 0.8GeV < pT

    // tpc signal
    histos.add("tpcSignal", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histos.add("tpcSignal_Pi", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

    // tof beta
    histos.add("tofBeta", "TOF beta vs p; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histos.add("tofBeta_Pi", "TOF beta vs p for pions; p [GeV/c]; #beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

    // Other signals
    histos.add("FT0A", "T0A amplitude", kTH1F, {{200, 0.0, 500.0}});
    histos.add("FT0C", "T0C amplitude", kTH1F, {{200, 0.0, 500.0}});
    histos.add("ZDC_A", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histos.add("ZDC_C", "ZDC amplitude", kTH1F, {{1000, 0.0, 15}});
    histos.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 100}});

    // Collin Soper Theta and Phi
    histos.add("CS_phi_pair_1", "#phi Distribution; #phi; Entries", kTH1F, {{nBins_Phi, -3.2, 3.2}});
    histos.add("CS_phi_pair_2", "#phi Distribution; #phi; Entries", kTH1F, {{nBins_Phi, -3.2, 3.2}});
    histos.add("CS_costheta_pair_1", "#theta Distribution;cos(#theta); Entries", kTH1F, {{nBins_CosTheta, -1, 1}});
    histos.add("CS_costheta_pair_2", "#theta Distribution;cos(#theta); Entries", kTH1F, {{nBins_CosTheta, -1, 1}});

  } // End of init function
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Calculate the Collins-Soper Frame----------------------------------------------------------------------------------------------------------------------------
  Double_t CosThetaCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    //  TVector3 beta = (-1. / fourpion.E()) * fourpion.Vect();
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(pTargCM).Vect()).Unit()};

    // Axes
    ROOT::Math::XYZVectorF zaxis_CS{((Beam1_CM.Unit() - Beam2_CM.Unit()).Unit())};

    Double_t CosThetaCS = zaxis_CS.Dot((v1_CM));
    return CosThetaCS;
  } // End of CosThetaCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  // Calculate Phi in Collins-Soper Frame------------------------------------------------------------------------------------------------------------------------
  Double_t PhiCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;
    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxis_CS{((Beam1_CM.Unit() - Beam2_CM.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
    ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};

    Double_t phi = TMath::ATan2(yaxis_CS.Dot(v1_CM), xaxis_CS.Dot(v1_CM));
    return phi;
  } // End of PhiCollinsSoperFrame function------------------------------------------------------------------------------------------------------------------------

  // Begin of Process function--------------------------------------------------------------------------------------------------------------------------------------------------
  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {

    int gapSide = collision.gapSide();
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    histos.fill(HIST("GapSide"), gapSide);
    histos.fill(HIST("TrueGapSide"), truegapSide);
    histos.fill(HIST("EventCounts"), 1);
    gapSide = truegapSide;

    if ((gapSide != 2)) {
      return;
    }

    histos.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histos.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histos.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histos.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histos.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    if (StrictEventSelection) {
      if (collision.numContrib() != 4) {
        return;
      }
    } else {
      if (collision.numContrib() >= 10) {
        return;
      }
    }

    std::vector<decltype(tracks.begin())> WOTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_tracks;
    std::vector<decltype(tracks.begin())> WTS_PID_Pi_tracks;
    std::vector<decltype(tracks.begin())> Pi_plus_tracks;
    std::vector<decltype(tracks.begin())> Pi_minus_tracks;

    for (auto& t0 : tracks) {

      WOTS_tracks.push_back(t0);

      if (trackselector(t0, parameters)) {
        WTS_tracks.push_back(t0);

        if (selectionPIDPion(t0, true, nSigmaTPC_cut, nSigmaTOF_cut)) {
          WTS_PID_Pi_tracks.push_back(t0);
          if (t0.sign() == 1) {
            Pi_plus_tracks.push_back(t0);
          }
          if (t0.sign() == -1) {
            Pi_minus_tracks.push_back(t0);
          }
        } // End of Selection PID Pion

      } // End of track selections

    } // End of loop over tracks

    int len_WOTS = static_cast<int>(WOTS_tracks.size());
    int len_WTS = static_cast<int>(WTS_tracks.size());
    int len_WTS_PID_Pi = static_cast<int>(WTS_PID_Pi_tracks.size());
    int len_Pi_plus = static_cast<int>(Pi_plus_tracks.size());
    int len_Pi_minus = static_cast<int>(Pi_minus_tracks.size());

    for (int i = 0; i < len_WOTS; i++) {
      histos.fill(HIST("tpcNSigmaPi_WOTS"), WOTS_tracks[i].tpcNSigmaPi());
      histos.fill(HIST("tofNSigmaPi_WOTS"), WOTS_tracks[i].tofNSigmaPi());
      histos.fill(HIST("pT_track_WOTS"), TMath::Sqrt(WOTS_tracks[i].px() * WOTS_tracks[i].px() + WOTS_tracks[i].py() * WOTS_tracks[i].py()));
    } // End of loop over tracks without selection

    for (int i = 0; i < len_WTS; i++) {

      histos.fill(HIST("tpcSignal"), TMath::Sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].tpcSignal());
      histos.fill(HIST("tofBeta"), TMath::Sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py() + WTS_tracks[i].pz() * WTS_tracks[i].pz()), WTS_tracks[i].beta());
      histos.fill(HIST("tpcNSigmaPi_WTS"), WTS_tracks[i].tpcNSigmaPi());
      histos.fill(HIST("tofNSigmaPi_WTS"), WTS_tracks[i].tofNSigmaPi());
      histos.fill(HIST("pT_track_WTS"), TMath::Sqrt(WTS_tracks[i].px() * WTS_tracks[i].px() + WTS_tracks[i].py() * WTS_tracks[i].py()));
    } // End of loop over tracks with selection only

    for (int i = 0; i < len_WTS_PID_Pi; i++) {

      histos.fill(HIST("tpcSignal_Pi"), TMath::Sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].tpcSignal());
      histos.fill(HIST("tofBeta_Pi"), TMath::Sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py() + WTS_PID_Pi_tracks[i].pz() * WTS_PID_Pi_tracks[i].pz()), WTS_PID_Pi_tracks[i].beta());

      histos.fill(HIST("tpcNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tpcNSigmaPi());
      histos.fill(HIST("tpcNSigmaPi_WTS_PID_Pi_Ka"), WTS_PID_Pi_tracks[i].tpcNSigmaKa());
      histos.fill(HIST("tpcNSigmaPi_WTS_PID_Pi_Pr"), WTS_PID_Pi_tracks[i].tpcNSigmaPr());
      histos.fill(HIST("tpcNSigmaPi_WTS_PID_Pi_El"), WTS_PID_Pi_tracks[i].tpcNSigmaEl());
      histos.fill(HIST("tpcNSigmaPi_WTS_PID_Pi_Mu"), WTS_PID_Pi_tracks[i].tpcNSigmaMu());

      histos.fill(HIST("tofNSigmaPi_WTS_PID_Pi"), WTS_PID_Pi_tracks[i].tofNSigmaPi());
      histos.fill(HIST("tofNSigmaPi_WTS_PID_Pi_Ka"), WTS_PID_Pi_tracks[i].tofNSigmaKa());
      histos.fill(HIST("tofNSigmaPi_WTS_PID_Pi_Pr"), WTS_PID_Pi_tracks[i].tofNSigmaPr());
      histos.fill(HIST("tofNSigmaPi_WTS_PID_Pi_El"), WTS_PID_Pi_tracks[i].tofNSigmaEl());
      histos.fill(HIST("tofNSigmaPi_WTS_PID_Pi_Mu"), WTS_PID_Pi_tracks[i].tofNSigmaMu());

      histos.fill(HIST("pT_track_WTS_PID_Pi"), TMath::Sqrt(WTS_PID_Pi_tracks[i].px() * WTS_PID_Pi_tracks[i].px() + WTS_PID_Pi_tracks[i].py() * WTS_PID_Pi_tracks[i].py()));
    } // End of loop over tracks with selection and PID selection of Pions

    if (len_WTS_PID_Pi != 4) {
      return;
    }

    // Selecting Events with net charge = 0
    if (len_Pi_minus == 2 && len_Pi_plus == 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      ROOT::Math::PtEtaPhiMVector k1, k2, k3, k4, k1234, k13, k14, k23, k24;

      p1.SetXYZM(Pi_plus_tracks[0].px(), Pi_plus_tracks[0].py(), Pi_plus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(Pi_plus_tracks[1].px(), Pi_plus_tracks[1].py(), Pi_plus_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(Pi_minus_tracks[0].px(), Pi_minus_tracks[0].py(), Pi_minus_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(Pi_minus_tracks[1].px(), Pi_minus_tracks[1].py(), Pi_minus_tracks[1].pz(), o2::constants::physics::MassPionCharged);

      k1.SetCoordinates(p1.Pt(), p1.Eta(), p1.Phi(), o2::constants::physics::MassPionCharged);
      k2.SetCoordinates(p2.Pt(), p2.Eta(), p2.Phi(), o2::constants::physics::MassPionCharged);
      k3.SetCoordinates(p3.Pt(), p3.Eta(), p3.Phi(), o2::constants::physics::MassPionCharged);
      k4.SetCoordinates(p4.Pt(), p4.Eta(), p4.Phi(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;
      k1234 = k1 + k2 + k3 + k4;

      k13 = k1 + k3;
      k14 = k1 + k4;
      k23 = k2 + k3;
      k24 = k2 + k4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histos.fill(HIST("pT_event_0charge_WTS_PID_Pi"), p1234.Pt());
        if (p1234.Pt() < 0.15) {
          histos.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainA"), p1234.M());

          auto phi_pair_1 = PhiCollinsSoperFrame(k13, k24, k1234);
          auto phi_pair_2 = PhiCollinsSoperFrame(k14, k23, k1234);
          auto cos_theta_1 = CosThetaCollinsSoperFrame(k13, k24, k1234);
          auto cos_theta_2 = CosThetaCollinsSoperFrame(k14, k23, k1234);

          histos.fill(HIST("CS_phi_pair_1"), phi_pair_1);
          histos.fill(HIST("CS_phi_pair_2"), phi_pair_2);
          histos.fill(HIST("CS_costheta_pair_1"), cos_theta_1);
          histos.fill(HIST("CS_costheta_pair_2"), cos_theta_2);
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histos.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histos.fill(HIST("rapidity_event_0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for 0 charge events

    // Selecting Events with net charge != 0 for estimation of background
    if (len_Pi_minus != 2 && len_Pi_plus != 2) {

      TLorentzVector p1, p2, p3, p4, p1234;
      p1.SetXYZM(WTS_PID_Pi_tracks[0].px(), WTS_PID_Pi_tracks[0].py(), WTS_PID_Pi_tracks[0].pz(), o2::constants::physics::MassPionCharged);
      p2.SetXYZM(WTS_PID_Pi_tracks[1].px(), WTS_PID_Pi_tracks[1].py(), WTS_PID_Pi_tracks[1].pz(), o2::constants::physics::MassPionCharged);
      p3.SetXYZM(WTS_PID_Pi_tracks[2].px(), WTS_PID_Pi_tracks[2].py(), WTS_PID_Pi_tracks[2].pz(), o2::constants::physics::MassPionCharged);
      p4.SetXYZM(WTS_PID_Pi_tracks[3].px(), WTS_PID_Pi_tracks[3].py(), WTS_PID_Pi_tracks[3].pz(), o2::constants::physics::MassPionCharged);

      p1234 = p1 + p2 + p3 + p4;

      if (std::fabs(p1234.Rapidity()) < 0.5) {
        histos.fill(HIST("pT_event_non0charge_WTS_PID_Pi"), p1234.Pt());

        if (p1234.Pt() < 0.15) {
          histos.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainA"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainA"), p1234.M());
        }
        if (p1234.Pt() > 0.15 && p1234.Pt() < 0.80) {
          histos.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainB"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainB"), p1234.M());
        }
        if (p1234.Pt() > 0.80) {
          histos.fill(HIST("rapidity_event_non0charge_WTS_PID_Pi_domainC"), p1234.Rapidity());
          histos.fill(HIST("invMass_event_non0charge_WTS_PID_Pi_domainC"), p1234.M());
        }
      } // End of Rapidity range selection

    } // End of Analysis for non 0 charge events

  } // End of process function
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // End of Struct UPCAnalysis
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<exclusiveRhoTo4Pi>(cfgc)};
}
