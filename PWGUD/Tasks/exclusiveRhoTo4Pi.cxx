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
struct UPCAnalysis {
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
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 85, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.5, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0, "Track Pt"};

  Configurable<float> nSigmaTPC_cut{"nsigmatpccut", 3, "TPC cut"};
  Configurable<float> nSigmaTOF_cut{"nsigmatofcut", 9, "TOF cut"};
  Configurable<bool> useTOF{"useTOF", true, "use TOF"};
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // Begin of Init Function-----------------------------------------------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {

    histos.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    histos.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    histos.add("EventsCounts", "Total Events; Entries", kTH1F, {{100, 0, 100}}); // 2=#Events, 3= #Selected Events, 5=#Selected Events with Zero Net charge, 6=Selected Events with Non-Zero Net charge
    histos.add("nSelectedEvents", "Number of Selected Events with 4 good tracks; Entries", kTH1F, {{10, 0, 10}});
    histos.add("nSelectedEventswithnonzerocharge", "Number of Selected Good Events with Non Zero Net charge; Entries", kTH1F, {{10, 0, 10}});
    histos.add("nSelectedEventswithzerocharge", "Number of Selected Good Events with Zero Net charge; Entries", kTH1F, {{10, 0, 10}});

    histos.add("tpc_nsigma_pion_WOTS", "TPC nSigma Pion; Entries", kTH1F, {{100, -15, 15}});
    histos.add("tpc_nsigma_electron_WOTS", "TPC nSigma Electron; Entries", kTH1F, {{100, -15, 15}});
    histos.add("tpc_nsigma_pion_electron_WOTS", "TPC nSigma Pion vs Electron; N sigma e-; N sigma Pi", kTH2F, {{300, -15, 15}, {300, -15, 15}});
    histos.add("tpc_nsigma_pion_eletron_WTS", "TPC nSigma Pion vs Electron with track selection; N sigma e-; N sigma Pi", kTH2F, {{300, -15, 15}, {300, -15, 15}});

    histos.add("tof_nsigma_pion_WOTS", "TOF nSigma Pion; Entries", kTH1F, {{100, -15, 15}});
    histos.add("tof_nsigma_electron_WOTS", "TOF nSigma Electron; Entries", kTH1F, {{100, -15, 15}});
    histos.add("tof_nsigma_pion_electron_WOTS", "TOF nSigma Pion vs Electron; N sigma e-; N sigma Pi", kTH2F, {{300, -15, 15}, {300, -15, 15}});
    histos.add("tof_nsigma_pion_eletron_WTS", "TOF nSigma Pion vs Electron with track selection; N sigma e-; N sigma Pi", kTH2F, {{300, -15, 15}, {300, -15, 15}});

    histos.add("pT_with_nsigma_3_WOTS", "pT with nsigma 3 Without track selection; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_with_nsigma_3_WTS", "pT with nsigma 3 with track selection; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_WTS", "pT with track selection; pT; Entries", kTH1F, {{100, 0, 2}});
    histos.add("pT_WOTS", "pT without track selection; pT; Entries", kTH1F, {{100, 0, 2}});

    histos.add("pT_event_0c_WOTS", "Event pT in 0 Charge Events Without Track Selection; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_event_n0c_WOTS", "Event pT in Non 0 Charge Events Without Track Selection; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});

    histos.add("pT_event_0c", "Event pT in 0 Charge Events; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_event_n0c", "Event pT in Non 0 Charge Events; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_pi_0c", "pT of Pions in 0 Charge Events; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});
    histos.add("pT_pi_n0c", "pT of Pions in Non 0 Charge Events; pT [GeV/c]; Counts", kTH1F, {{100, 0, 2}});

    // invariant mass of 0 charge events
    histos.add("invmass_0c_domA", "Invariant Mass Distribution of 0 charge Events in the region pT<0.15 GeV/c", kTH1F, {{1000, 0.8, 2.5}});
    histos.add("invmass_0c_domB", "Invariant Mass Distribution of 0 charge Events in the region 0.15 < pT < 0.8 GeV/c", kTH1F, {{100, 0.8, 2.5}});
    histos.add("invmass_0c_domC", "Invariant Mass Distribution of 0 charge Events in the region pT>0.8 GeV/c", kTH1F, {{100, 0.8, 2.5}});

    // invariant mass of non 0 charge events
    histos.add("invmass_n0c_domA", "Invariant Mass Distribution of non 0 charge events in the region pT<0.15 GeV/c", kTH1F, {{100, 0.8, 2.5}});
    histos.add("invmass_n0c_domB", "Invariant Mass Distribution of non 0 charge events in the region 0.15 < pT < 0.8 GeV/c", kTH1F, {{100, 0.8, 2.5}});
    histos.add("invmass_n0c_domC", "Invariant Mass Distribution of non 0 charge events in the region pT>0.8 GeV/c", kTH1F, {{100, 0.8, 2.5}});

    // invariant mass of like sign and unlike sign pairs
    histos.add("invmass_unlikesign_event", "Event Mass Distribution of unlike sign pairs", kTH2F, {{100, 0.8, 2.5}, {100, 0.8, 2.5}});
    histos.add("invmass_likesign_event", "Event Mass Distribution of like sign pairs", kTH2F, {{100, 0.8, 2.5}, {100, 0.8, 2.5}});
    histos.add("invmass_total", "Invariant Mass Distribution of the Event", kTH1F, {{100, 0.8, 2.5}});

    // QA plots

    // tpc signal
    histos.add("tpc_dEdx", "TPC dEdx vs p; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histos.add("tpc_dEdx_pion", "TPC dEdx vs p for pions; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histos.add("tpc_dEdx_kaons", "TPC dEdx vs p for ekaons; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});
    histos.add("tpc_dEdx_protons", "TPC dEdx vs p for protons; p [GeV/c]; dEdx [a.u.]", kTH2F, {{500, 0, 10}, {5000, 0.0, 5000.0}});

    // tof beta
    histos.add("tof_beta", "TOF beta vs p; p [GeV/c]; beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histos.add("tof_beta_pion", "TOF beta vs p for pions; p [GeV/c]; beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histos.add("tof_beta_kaons", "TOF beta vs p for kaons; p [GeV/c]; beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});
    histos.add("tof_beta_protons", "TOF beta vs p for protons; p [GeV/c]; beta", kTH2F, {{500, 0, 10}, {500, 0.0, 1.0}});

    // Other signals
    histos.add("FT0A", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
    histos.add("FT0C", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
    histos.add("ZDC_A", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
    histos.add("ZDC_C", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
    histos.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});

    // Collin Sopen Stuff
    histos.add("CSphi", "#phi Distribution; #phi; Entries", kTH1F, {{200, -3.2, 3.2}});
    histos.add("CStheta", "#theta Distribution; #theta; Entries", kTH1F, {{200, -3.2, 3.2}});
    histos.add("CS_Theta_vs_phi", "#theta vs #phi; #phi; #theta", kTH2F, {{200, -3.2, 3.2}, {200, 0.0, 3.2}});

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
    histos.fill(HIST("EventsCounts"), 1);
    gapSide = truegapSide;

    if ((gapSide != 2)) {
      return;
    }

    if (collision.numContrib() != 4) {
      return;
    }

    histos.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
    histos.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
    histos.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
    histos.fill(HIST("ZDC_A"), collision.energyCommonZNA());
    histos.fill(HIST("ZDC_C"), collision.energyCommonZNC());

    std::vector<decltype(tracks.begin())> selectedTracks_p;          // positive charge tracks
    std::vector<decltype(tracks.begin())> selectedTracks_n;          // negative charge tracks
    std::vector<decltype(tracks.begin())> allSelectedTracks;         // all selected tracks
    std::vector<decltype(tracks.begin())> allTrackswithoutselection; // all tracks without selection

    for (auto& t0 : tracks) {

      // TPC signal without track selection
      histos.fill(HIST("tpc_nsigma_pion_WOTS"), t0.tpcNSigmaPi());
      histos.fill(HIST("tpc_nsigma_electron_WOTS"), t0.tpcNSigmaEl());
      histos.fill(HIST("tpc_nsigma_pion_electron_WOTS"), t0.tpcNSigmaEl(), t0.tpcNSigmaPi());

      // TOF signal without track selection
      histos.fill(HIST("tof_nsigma_pion_WOTS"), t0.tofNSigmaPi());
      histos.fill(HIST("tof_nsigma_electron_WOTS"), t0.tofNSigmaEl());
      histos.fill(HIST("tof_nsigma_pion_electron_WOTS"), t0.tofNSigmaEl(), t0.tofNSigmaPi());

      histos.fill(HIST("pT_WOTS"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py()));

      // pT with nsigma 3 without track selection
      if (fabs(t0.tpcNSigmaPi()) < 3.0 && fabs(t0.tpcNSigmaEl()) < 3.0) {
        histos.fill(HIST("pT_with_nsigma_3_WOTS"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py()));
      }

      allTrackswithoutselection.push_back(t0);
      // track selection
      if (!trackselector(t0, parameters)) {
        continue;
      }

      // dE/dx for selected tracks

      histos.fill(HIST("tpc_dEdx"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.tpcSignal());
      histos.fill(HIST("tof_beta"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.beta());

      // pT with nsigma 3 with track selection
      if (fabs(t0.tpcNSigmaPi()) < 3.0 && fabs(t0.tpcNSigmaEl()) < 3.0) {
        histos.fill(HIST("pT_with_nsigma_3_WTS"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py()));
      }

      // pT with track selection
      histos.fill(HIST("pT_WTS"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py()));

      if ((selectionPIDKaon(t0, true, nSigmaTPC_cut, nSigmaTOF_cut))) {
        histos.fill(HIST("tpc_dEdx_kaons"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.tpcSignal());
        histos.fill(HIST("tof_beta_kaons"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.beta());
      }

      if ((selectionPIDProton(t0, true, nSigmaTPC_cut, nSigmaTOF_cut))) {
        histos.fill(HIST("tpc_dEdx_protons"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.tpcSignal());
        histos.fill(HIST("tof_beta_protons"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.beta());
      }

      // TPC signal with track selection
      histos.fill(HIST("tpc_nsigma_pion_eletron_WTS"), t0.tpcNSigmaEl(), t0.tpcNSigmaPi());
      histos.fill(HIST("tof_nsigma_pion_eletron_WTS"), t0.tofNSigmaEl(), t0.tofNSigmaPi());

      if (!(selectionPIDPion(t0, true, nSigmaTPC_cut, nSigmaTOF_cut))) {
        continue;
      }
      histos.fill(HIST("tpc_dEdx_pion"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.tpcSignal());
      histos.fill(HIST("tof_beta_pion"), TMath::Sqrt(t0.px() * t0.px() + t0.py() * t0.py() + t0.pz() * t0.pz()), t0.beta());

      allSelectedTracks.push_back(t0);
      if (t0.sign() > 0) {
        selectedTracks_p.push_back(t0);
      } else if (t0.sign() < 0) {
        selectedTracks_n.push_back(t0);
      }
    } // End of loop over tracks

    if (allTrackswithoutselection.size() == 4) {

      histos.fill(HIST("nSelectedEvents"), 4);
      TLorentzVector v1, v2, v3, v4, v1234;
      v1.SetXYZM(allTrackswithoutselection[0].px(), allTrackswithoutselection[0].py(), allTrackswithoutselection[0].pz(), o2::constants::physics::MassPionCharged);
      v2.SetXYZM(allTrackswithoutselection[1].px(), allTrackswithoutselection[1].py(), allTrackswithoutselection[1].pz(), o2::constants::physics::MassPionCharged);
      v3.SetXYZM(allTrackswithoutselection[2].px(), allTrackswithoutselection[2].py(), allTrackswithoutselection[2].pz(), o2::constants::physics::MassPionCharged);
      v4.SetXYZM(allTrackswithoutselection[3].px(), allTrackswithoutselection[3].py(), allTrackswithoutselection[3].pz(), o2::constants::physics::MassPionCharged);
      int charge = allTrackswithoutselection[0].sign() + allTrackswithoutselection[1].sign() + allTrackswithoutselection[2].sign() + allTrackswithoutselection[3].sign();
      v1234 = v1 + v2 + v3 + v4;
      if (charge == 0) {
        histos.fill(HIST("pT_event_0c_WOTS"), v1234.Pt());
      } else {
        histos.fill(HIST("pT_event_n0c_WOTS"), v1234.Pt());
      }
    } // end of if for 4 track events without selection

    static_cast<int> totalSize = int(selectedTracks_p.size() + selectedTracks_n.size());

    if ((selectedTracks_p.size() == 2) && (selectedTracks_n.size() == 2)) {

      histos.fill(HIST("EventsCounts"), 3);
      histos.fill(HIST("EventsCounts"), 5);

      TLorentzVector v1, v2, v3, v4, v1234, v12, v34, v13, v24, v14, v23;
      v1.SetXYZM(selectedTracks_p[0].px(), selectedTracks_p[0].py(), selectedTracks_p[0].pz(), o2::constants::physics::MassPionCharged); // Pi+
      v2.SetXYZM(selectedTracks_p[1].px(), selectedTracks_p[1].py(), selectedTracks_p[1].pz(), o2::constants::physics::MassPionCharged); // Pi+

      v3.SetXYZM(selectedTracks_n[0].px(), selectedTracks_n[0].py(), selectedTracks_n[0].pz(), o2::constants::physics::MassPionCharged); // Pi-
      v4.SetXYZM(selectedTracks_n[1].px(), selectedTracks_n[1].py(), selectedTracks_n[1].pz(), o2::constants::physics::MassPionCharged); // Pi-

      v1234 = v1 + v2 + v3 + v4;

      v12 = v1 + v2;
      v34 = v3 + v4;
      v13 = v1 + v3;
      v24 = v2 + v4;
      v14 = v1 + v4;
      v23 = v2 + v3;

      ROOT::Math::PtEtaPhiMVector v1234_1(v1234.Pt(), v1234.Eta(), v1234.Phi(), o2::constants::physics::MassPionCharged);

      ROOT::Math::PtEtaPhiMVector v13_1(v13.Pt(), v13.Eta(), v13.Phi(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PtEtaPhiMVector v24_1(v24.Pt(), v24.Eta(), v24.Phi(), o2::constants::physics::MassPionCharged);

      ROOT::Math::PtEtaPhiMVector v14_1(v14.Pt(), v14.Eta(), v14.Phi(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PtEtaPhiMVector v23_1(v23.Pt(), v23.Eta(), v23.Phi(), o2::constants::physics::MassPionCharged);

      histos.fill(HIST("CSphi"), PhiCollinsSoperFrame(v13_1, v24_1, v1234_1));
      histos.fill(HIST("CStheta"), CosThetaCollinsSoperFrame(v13_1, v24_1, v1234_1));
      histos.fill(HIST("CS_Theta_vs_phi"), PhiCollinsSoperFrame(v13_1, v24_1, v1234_1), CosThetaCollinsSoperFrame(v13_1, v24_1, v1234_1));

      histos.fill(HIST("CSphi"), PhiCollinsSoperFrame(v14_1, v23_1, v1234_1));
      histos.fill(HIST("CStheta"), CosThetaCollinsSoperFrame(v14_1, v23_1, v1234_1));
      histos.fill(HIST("CS_Theta_vs_phi"), PhiCollinsSoperFrame(v14_1, v23_1, v1234_1), CosThetaCollinsSoperFrame(v14_1, v23_1, v1234_1));

      if (fabs(v1234.Rapidity()) < 0.5) {
        histos.fill(HIST("nSelectedEventswithzerocharge"), 6);
        histos.fill(HIST("pT_pi_0c"), v1.Pt());
        histos.fill(HIST("pT_pi_0c"), v2.Pt());
        histos.fill(HIST("pT_pi_0c"), v3.Pt());
        histos.fill(HIST("pT_pi_0c"), v4.Pt());
        histos.fill(HIST("pT_event_0c"), v1234.Pt());

        if (v1234.Pt() < 0.15) {
          histos.fill(HIST("invmass_0c_domA"), v1234.M());
        } else if (v1234.Pt() > 0.15 && v1234.Pt() < 0.8) {
          histos.fill(HIST("invmass_0c_domB"), v1234.M());
        } else if (v1234.Pt() > 0.8) {
          histos.fill(HIST("invmass_0c_domC"), v1234.M());
        }

        histos.fill(HIST("invmass_unlikesign_event"), v13.M(), v24.M());
        histos.fill(HIST("invmass_unlikesign_event"), v14.M(), v23.M());
        histos.fill(HIST("invmass_likesign_event"), v12.M(), v34.M());
        if (v1234.Pt() < 2.0) {
          histos.fill(HIST("invmass_total"), v1234.M());
        }

      } // end of rapidity cut

    } // end of if for 2 positive and 2 negative charge events

    if (totalSize == 4 && static_cast<int>(allSelectedTracks.size()) == 4) {

      histos.fill(HIST("EventsCounts"), 3);
      histos.fill(HIST("EventsCounts"), 7);

      histos.fill(HIST("nSelectedEventswithnonzerocharge"), 4);
      TLorentzVector v1, v2, v3, v4, v1234;
      v1.SetXYZM(allSelectedTracks[0].px(), allSelectedTracks[0].py(), allSelectedTracks[0].pz(), o2::constants::physics::MassPionCharged);
      v2.SetXYZM(allSelectedTracks[1].px(), allSelectedTracks[1].py(), allSelectedTracks[1].pz(), o2::constants::physics::MassPionCharged);
      v3.SetXYZM(allSelectedTracks[2].px(), allSelectedTracks[2].py(), allSelectedTracks[2].pz(), o2::constants::physics::MassPionCharged);
      v4.SetXYZM(allSelectedTracks[3].px(), allSelectedTracks[3].py(), allSelectedTracks[3].pz(), o2::constants::physics::MassPionCharged);

      v1234 = v1 + v2 + v3 + v4;

      if (fabs(v1234.Rapidity()) < 0.5) {
        histos.fill(HIST("pT_pi_n0c"), v1.Pt());
        histos.fill(HIST("pT_pi_n0c"), v2.Pt());
        histos.fill(HIST("pT_pi_n0c"), v3.Pt());
        histos.fill(HIST("pT_pi_n0c"), v4.Pt());
        histos.fill(HIST("pT_event_n0c"), v1234.Pt());

        if (v1234.Pt() < 0.15) {
          histos.fill(HIST("invmass_n0c_domA"), v1234.M());
        } else if (v1234.Pt() > 0.15 && v1234.Pt() < 0.8) {
          histos.fill(HIST("invmass_n0c_domB"), v1234.M());
        } else if (v1234.Pt() > 0.8) {
          histos.fill(HIST("invmass_n0c_domC"), v1234.M());
        }
      } // end of rapidity cut

    } // end else if for non zero charge events

  } // End of process function
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // End of Struct UPCAnalysis
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCAnalysis>(cfgc)};
}
