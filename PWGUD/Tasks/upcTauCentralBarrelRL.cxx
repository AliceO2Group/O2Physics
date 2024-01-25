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
/// \brief
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  12.07.2022

// #include <algorithm>
// #include <iterator>

// O2 headers
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TEfficiency.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcTauCentralBarrelRL {

  // Global varialbes
  bool isFirstReconstructedCollisions;
  int countCollisions;
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry histos{
    "histos",
    {{"Events/hCountCollisions", ";Number of analysed collision (-)", {HistType::kTH1D, {{10, 0.5, 10.5}}}},
     {"Events/hNreconstructedTracks", ";Number of tracks in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedPVGTelectrons", ";Number of good track electrons from primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedPVGTmuons", ";Number of good track muons from primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedPVGTpions", ";Number of good track pions from primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedPVGTothers", ";Number of good track NOT electron/muon/pion particles from primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedPVGT", ";Number of good track particles from primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hNreconstructedNotPVGT", ";Number of good track particles from NOT primary vertex in a collision (-);Number of events (-)", {HistType::kTH1D, {{30, -0.5, 29.5}}}},
     {"Events/hChannelsRatio", ";Channels (-);Branching Ratio (-)", {HistType::kTH1D, {{10, -0.5, 9.5}}}}}};
  HistogramRegistry histosPID{
    "histosPID",
    {{"Tracks/raw/PID/hTPCsignalVsZ", "All tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/hTPCsignalVsP", "All tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/hTPCsignalVsPt", "All tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/hTPCsignalVsEta", "All tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/hTPCsignalVsPhi", "All tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/raw/PID/PosCharge/hTPCsignalVsZ", "Positively charged track;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/PosCharge/hTPCsignalVsP", "Positively charged track;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/PosCharge/hTPCsignalVsPt", "Positively charged track;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/PosCharge/hTPCsignalVsEta", "Positively charged track;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/PosCharge/hTPCsignalVsPhi", "Positively charged track;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/raw/PID/NegCharge/hTPCsignalVsZ", "Negatively charged track;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/NegCharge/hTPCsignalVsP", "Negatively charged track;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/NegCharge/hTPCsignalVsPt", "Negatively charged track;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/NegCharge/hTPCsignalVsEta", "Negatively charged track;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/raw/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged track;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/hTPCsignalVsZ", "All good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/hTPCsignalVsP", "All good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/hTPCsignalVsPt", "All good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/hTPCsignalVsEta", "All good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/hTPCsignalVsPhi", "All good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsZ", "Positively charged track;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP", "Positively charged track;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt", "Positively charged track;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta", "Positively charged track;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi", "Positively charged track;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsZ", "Negatively charged track;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP", "Negatively charged track;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt", "Negatively charged track;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta", "Negatively charged track;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged track;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Electron/hTPCsignalVsZ", "Identified electron;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Electron/hTPCsignalVsP", "Identified electron;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt", "Identified electron;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta", "Identified electron;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi", "Identified electron;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Muon/hTPCsignalVsZ", "Identified Muon;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Muon/hTPCsignalVsP", "Identified Muon;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt", "Identified Muon;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta", "Identified Muon;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi", "Identified Muon;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Pion/hTPCsignalVsZ", "Identified Pion;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Pion/hTPCsignalVsP", "Identified Pion;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt", "Identified Pion;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta", "Identified Pion;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi", "Identified Pion;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Others/hTPCsignalVsZ", "Identified NOT electron/Muon/Pion;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, -20., 20.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Others/hTPCsignalVsP", "Identified NOT electron/Muon/Pion;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Others/hTPCsignalVsPt", "Identified NOT electron/Muon/Pion;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Others/hTPCsignalVsEta", "Identified NOT electron/Muon/Pion;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{500, -2., 2.}, {200, 0., 200}}}},
     {"Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi", "Identified NOT electron/Muon/Pion;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{64, 0, 2 * TMath::Pi()}, {200, 0., 200}}}},
     {"EventTwoTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/TwoMuons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/TwoPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/ElectronPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventTwoTracks/MuonPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventFourTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventFourTracks/WithElectron/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventFourTracks/WithMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventFourTracks/WithPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventSixTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}},
     {"EventSixTracks/SixPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", {HistType::kTH2D, {{200, 0., 2.}, {200, 0., 200}}}}}};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", true, {"Print general info to terminal; default it true."}};
  Configurable<bool> verboseDebug{"verboseDebug", false, {"Print debug info to terminal; default it false."}};
  Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};

  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
  using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;
  using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions>::iterator;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printLargeMessage("INIT METHOD");
    countCollisions = 0;
    isFirstReconstructedCollisions = true;

    const AxisSpec axisZvtx{40, -20., 20.};
    const AxisSpec axisInvMass{400, 1., 5.};
    const AxisSpec axisInvMassWide{500, 0., 10.};
    const AxisSpec axisMom{40, 0., 2.};
    const AxisSpec axisMomWide{100, 0., 10.};
    const AxisSpec axisPt{40, 0., 2.};
    const AxisSpec axisPhi{64, -2 * TMath::Pi(), 2 * TMath::Pi()};
    const AxisSpec axisEta{50, -1.2, 1.2};
    const AxisSpec axisRap{50, -1.2, 1.2};
    const AxisSpec axisAcoplanarity{32, 0.0, o2::constants::math::PI};

    histos.add("Tracks/raw/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {axisZvtx});
    histos.add("Tracks/raw/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("Tracks/raw/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("Tracks/raw/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("Tracks/raw/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {axisEta});

    histos.add("Tracks/GoodTrack/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {axisZvtx});
    histos.add("Tracks/GoodTrack/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("Tracks/GoodTrack/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("Tracks/GoodTrack/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("Tracks/GoodTrack/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {axisEta});

    histos.add("EventTwoTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/hInvariantMassWideAllPionMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/hInvariantMassWideAllPionMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/TwoElectrons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/TwoElectrons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/TwoElectrons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/TwoElectrons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/TwoElectrons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/TwoElectrons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/TwoElectrons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/TwoElectrons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/TwoElectrons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/TwoElectrons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/TwoElectrons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/TwoElectrons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/TwoElectrons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
    histos.add("EventTwoTracks/TwoElectrons/hLeadingP", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/TwoElectrons/hLeadingPwide", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/TwoElectrons/hLeadingPt", ";Leading #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/TwoElectrons/hLeadingPhi", ";Leading #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/TwoElectrons/hLeadingRapidity", ";Leading #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventTwoTracks/TwoMuons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/TwoMuons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/TwoMuons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/TwoMuons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/TwoMuons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/TwoMuons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/TwoMuons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/TwoMuons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/TwoMuons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/TwoMuons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/TwoMuons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/TwoMuons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/TwoMuons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/TwoPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/TwoPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/TwoPions/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/TwoPions/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/TwoPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/TwoPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/TwoPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/TwoPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/TwoPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/TwoPions/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/TwoPions/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/TwoPions/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/TwoPions/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/TwoPions/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/ElectronMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/ElectronMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/ElectronMuon/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/ElectronMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/ElectronMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/ElectronMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/ElectronMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/ElectronMuon/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/ElectronMuon/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/ElectronMuon/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/ElectronMuon/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/ElectronMuon/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/ElectronPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/ElectronPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/ElectronPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/ElectronPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/ElectronPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/ElectronPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/ElectronPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/ElectronPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/ElectronPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/ElectronPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/ElectronPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/ElectronPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/MuonPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/MuonPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/MuonPion/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/MuonPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/MuonPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/MuonPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/MuonPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/MuonPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/MuonPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/MuonPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/MuonPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/MuonPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/MuonPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/MuonPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/ElectronMuPi/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
    histos.add("EventTwoTracks/ElectronMuPi/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/ElectronMuPi/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/ElectronMuPi/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/ElectronMuPi/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/ElectronMuPi/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronMuPi/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/ElectronMuPi/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/ElectronMuPi/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/ElectronMuPi/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronMuPi/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/ElectronMuPi/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/ElectronOther/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
    histos.add("EventTwoTracks/ElectronOther/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/ElectronOther/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/ElectronOther/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/ElectronOther/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/ElectronOther/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronOther/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/ElectronOther/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/ElectronOther/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/ElectronOther/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/ElectronOther/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/ElectronOther/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/ElectronOther/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/ElectronOther/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/ElectronOther/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventTwoTracks/PionsSelection/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/PionsSelection/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/PionsSelection/hInvariantMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventTwoTracks/PionsSelection/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
    histos.add("EventTwoTracks/PionsSelection/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventTwoTracks/PionsSelection/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventTwoTracks/PionsSelection/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventTwoTracks/PionsSelection/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventTwoTracks/PionsSelection/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
    histos.add("EventTwoTracks/PionsSelection/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
    histos.add("EventTwoTracks/PionsSelection/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
    histos.add("EventTwoTracks/PionsSelection/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
    histos.add("EventTwoTracks/PionsSelection/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add("EventTwoTracks/PionsSelection/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});

    histos.add("EventFourTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventFourTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventFourTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventFourTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventFourTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventFourTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventFourTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventFourTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventFourTracks/WithElectron/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventFourTracks/WithElectron/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventFourTracks/WithElectron/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventFourTracks/WithElectron/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventFourTracks/WithElectron/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventFourTracks/WithElectron/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventFourTracks/WithElectron/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventFourTracks/WithMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventFourTracks/WithMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventFourTracks/WithMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventFourTracks/WithMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventFourTracks/WithMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventFourTracks/WithMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventFourTracks/WithMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventFourTracks/WithPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventFourTracks/WithPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventFourTracks/WithPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventFourTracks/WithPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventFourTracks/WithPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventFourTracks/WithPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventFourTracks/WithPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventSixTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventSixTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventSixTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventSixTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventSixTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventSixTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventSixTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventSixTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

    histos.add("EventSixTracks/SixPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
    histos.add("EventSixTracks/SixPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
    histos.add("EventSixTracks/SixPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("EventSixTracks/SixPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
    histos.add("EventSixTracks/SixPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("EventSixTracks/SixPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("EventSixTracks/SixPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});

  } // end init

  // run (always called before process :( )
  void run(ProcessingContext& context)
  {

    if (verboseInfo)
      printLargeMessage("RUN METHOD");
    if (verboseDebug)
      LOGF(info, "countCollisions = %d", countCollisions);

  } // end run

  // process

  template <typename C, typename Ts>
  void fillHistograms(C reconstructedCollision, Ts reconstructedBarrelTracks)
  {

    if (isFirstReconstructedCollisions) {
      isFirstReconstructedCollisions = false;
      if (verboseInfo)
        printLargeMessage("START LOOPING OVER RECONSTRUCTED COLLISIONS");
    }

    histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(4); // 1, 2 and 3 are reserved for single-gap

    // Loop over tracks without selections
    for (auto& track : reconstructedBarrelTracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/raw/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/raw/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/raw/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/raw/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH1>(HIST("Tracks/raw/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      //          histosPID.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
      if (track.hasTPC()) {
        if (track.sign() == 1) {
          //                  histosPID.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else if (track.sign() == -1) {
          //                  histosPID.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
    } // Loop over tracks without selections

    int countPVGT = 0;
    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTothers = 0;
    int countPVGTpionsSelection = 0;
    std::vector<int> vecPVidx;
    // Loop over tracks with selections
    for (auto& track : reconstructedBarrelTracks) {
      if (track.isPVContributor() != 1)
        continue;
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      //          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
      histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
      if (track.hasTPC()) {
        if (track.sign() == 1) {
          //                  histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else if (track.sign() == -1) {
          //                  histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
      countPVGT++;
      int hypothesisID = testPIDhypothesis(track);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
          //                  histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
          //                  histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        } else {
          countPVGTpions++;
          //                  histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        }
      } else {
        countPVGTothers++;
        //              histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
        histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histosPID.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
      }
      if (abs(track.tpcNSigmaPi()) < 3)
        countPVGTpionsSelection++;

    } // Loop over tracks with selections

    histos.get<TH1>(HIST("Events/hNreconstructedPVGT"))->Fill(countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedNotPVGT"))->Fill(reconstructedBarrelTracks.size() - countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTelectrons"))->Fill(countPVGTelectrons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTmuons"))->Fill(countPVGTmuons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTpions"))->Fill(countPVGTpions);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTothers"))->Fill(countPVGTothers);

    if (countPVGTselected == 2) {
      TLorentzVector mother, daug[2], motherOfPions, pion[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      mother = daug[0] + daug[1];
      pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(211), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(211), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      motherOfPions = pion[0] + pion[1];
      auto acoplanarity = calculateAcoplanarity(daug[0].Phi(), daug[1].Phi());
      auto sign = trkDaug1.sign() * trkDaug2.sign();

      if (trkDaug1.hasTPC()) {
        histosPID.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTmuons == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histosPID.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTmuons == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 2)
          histosPID.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }

      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMass"))->Fill(motherOfPions.M());
      histos.get<TH1>(HIST("EventTwoTracks/hAcoplanarity"))->Fill(acoplanarity);
      histos.get<TH1>(HIST("EventTwoTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      if (motherOfPions.Pt() < 0.2) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassPtCut"))->Fill(motherOfPions.M());
      }

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTelectrons == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(0);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingRapidity"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()));
      }
      if (countPVGTmuons == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(1);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut"))->Fill(mother.M());
        }
      }
      if (countPVGTelectrons == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(2);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(3);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWidePtCut"))->Fill(mother.M());
        }
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(4);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(5);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWidePtCut"))->Fill(mother.M());
        }
      }
      if ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
        double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(12);
      }
      if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
        double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        if (countPVGTelectrons == 2)
          electronPt = (daug[0].Pt() > daug[1].Pt()) ? daug[0].Pt() : daug[1].Pt();
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(12);
      }
      if (countPVGTpionsSelection == 2) {
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMass"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWide"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherP"))->Fill(motherOfPions.P());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPwide"))->Fill(motherOfPions.P());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPt"))->Fill(motherOfPions.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPhi"))->Fill(motherOfPions.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherRapidity"))->Fill(motherOfPions.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersP"))->Fill(pion[0].P(), pion[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPwide"))->Fill(pion[0].P(), pion[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPt"))->Fill(pion[0].Pt(), pion[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPhi"))->Fill(pion[0].Phi(), pion[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersRapidity"))->Fill(pion[0].Rapidity(), pion[1].Rapidity());
        if (motherOfPions.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassPtCut"))->Fill(motherOfPions.M());
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCut"))->Fill(motherOfPions.M());
          if (sign < 0)
            histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUS"))->Fill(motherOfPions.M());
          if (sign > 0)
            histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLS"))->Fill(motherOfPions.M());
        }
      }
    } else if (countPVGTselected == 4) {

      TLorentzVector mother, daug[4];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3];

      if (trkDaug1.hasTPC()) {
        histosPID.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 4)
          histosPID.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histosPID.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histosPID.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 4)
          histosPID.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histosPID.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histosPID.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 4)
          histosPID.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histosPID.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histosPID.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 4)
          histosPID.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histosPID.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histosPID.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }

      histos.get<TH1>(HIST("EventFourTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventFourTracks/hMotherRapidity"))->Fill(mother.Rapidity());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 4) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(6);
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherRapidity"))->Fill(mother.Rapidity());
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 3) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(7);
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherRapidity"))->Fill(mother.Rapidity());
      }
      if (countPVGTpions == 3 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(8);
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherRapidity"))->Fill(mother.Rapidity());
      }
    } else if (countPVGTselected == 6) {
      TLorentzVector mother, daug[6];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      const auto& trkDaug5 = reconstructedBarrelTracks.iteratorAt(vecPVidx[4]);
      const auto& trkDaug6 = reconstructedBarrelTracks.iteratorAt(vecPVidx[5]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      daug[4].SetPxPyPzE(trkDaug5.px(), trkDaug5.py(), trkDaug5.pz(), energy(pdg->Mass(trackPDG(trkDaug5)), trkDaug5.px(), trkDaug5.py(), trkDaug5.pz()));
      daug[5].SetPxPyPzE(trkDaug6.px(), trkDaug6.py(), trkDaug6.pz(), energy(pdg->Mass(trackPDG(trkDaug6)), trkDaug6.px(), trkDaug6.py(), trkDaug6.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3] + daug[4] + daug[5];

      if (trkDaug1.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }
      if (trkDaug5.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
      }
      if (trkDaug6.hasTPC()) {
        histosPID.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
        if (countPVGTpions == 6)
          histosPID.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
      }

      histos.get<TH1>(HIST("EventSixTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventSixTracks/hMotherRapidity"))->Fill(mother.Rapidity());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 6) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(9);
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherRapidity"))->Fill(mother.Rapidity());
      }
    } else {
      if (verboseDebug) {
        printMediumMessage("Other particles");
      }
    }

  } // end fillHistograms

  void processDGrecoLevel(FullUDCollision const& reconstructedCollision,
                          FullUDTracks const& reconstructedBarrelTracks)
  {
    countCollisions++;

    fillHistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDGrecoLevel

  void processSGrecoLevel(FullSGUDCollision const& reconstructedCollision,
                          FullUDTracks const& reconstructedBarrelTracks)
  {
    countCollisions++;

    if (reconstructedCollision.gapSide() == 0)
      histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(1);
    if (reconstructedCollision.gapSide() == 1)
      histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(2);
    if (reconstructedCollision.gapSide() == 2)
      histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(3);
    if (reconstructedCollision.gapSide() != whichGapSide)
      return;

    fillHistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDGrecoLevel

  void processAnalysisFinished(FullUDCollision const& collisions)
  {

    if (verboseInfo)
      LOGF(info, "####################################### END OF THIS DATAFRAME #######################################");
    if (verboseDebug)
      LOGF(info, "countCollisions = %d");
    isFirstReconstructedCollisions = true;

  } // end processAnalysisFinished

  PROCESS_SWITCH(UpcTauCentralBarrelRL, processDGrecoLevel, "Iterate UD tables with reconstructed data created by DG-Candidate-Producer", false);
  PROCESS_SWITCH(UpcTauCentralBarrelRL, processSGrecoLevel, "Iterate UD tables with reconstructed data created by SG-Candidate-Producer", false);
  PROCESS_SWITCH(UpcTauCentralBarrelRL, processAnalysisFinished, "Simply runs in the end of the dataframe", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTauCentralBarrelRL>(cfgc, TaskName{"upc-tau-rl"})};
}
