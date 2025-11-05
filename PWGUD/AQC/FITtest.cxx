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
/// \brief A task for testing FIT selection for Ultra-perimpheral Collisions
/// \author Anisa Khatun, anisa.khatun@cern.ch
/// \since  04.08.2023

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/BCRange.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct FITtest {

  // inivinitialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::TOFSignal, aod::pidTOFbeta>;
  using FWs = aod::FwdTracks;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMain")) {

      // collisions
      registry.add("All/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
      registry.add("All/AllRelBC", "Total Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("All/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("All/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("All/PVCFIT", "PV contributors with FIT; PV contributors; Tracks", {HistType::kTH1F, {{100, -0.5, 99.5}}});
      registry.add("All/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("All/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("All/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("All/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
      registry.add("All/trketa", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("All/trkpt", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});

      // FV0
      registry.add("All/FV0/hV0A", "Time FV0A; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
      registry.add("All/FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/FV0/FV0A", "#FV0A; Channel; FV0A Amplitude", {HistType::kTH2F, {{48, -0.5, 47.5}, {2000, 0., 2000.}}});

      // FT0
      registry.add("All/FT0/hT0AC", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("All/FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -4.5, 1995.5}}});
      registry.add("All/FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -4.5, 1995.5}}});
      registry.add("All/FT0/FT0A", "#FT0A; Channel; FT0A Amplitude", {HistType::kTH2F, {{96, -0.5, 95.5}, {1000, 0., 1000.}}});
      registry.add("All/FT0/FT0C", "#FT0C; Channel; FT0C Amplitude", {HistType::kTH2F, {{112, -0.5, 111.5}, {1000, 0., 1000.}}});
      registry.add("All/FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // FDD
      registry.add("All/FDD/hFDDAC", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("All/FDD/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("All/FDD/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("All/FDD/FDDA", "#FDDA; Channel; FDDA Amplitude", {HistType::kTH2F, {{8, -0.5, 7.5}, {1000, 0., 1000.}}});
      registry.add("All/FDD/FDDC", "#FDDC; Channel; FDDC Amplitude", {HistType::kTH2F, {{8, -0.5, 7.5}, {1000, 0., 1000.}}});
      registry.add("All/FDD/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // FIT
      registry.add("All/FITAamp", "#FIT A side; FITA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("All/FITCamp", "#FIT C side; FITC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("All/FITACCorr", "FIT amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // ZDC

      registry.add("All/ZDC/hZNAC", "Time Correlation ZN; ZNA Time (ns) ; ZNC Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("All/ZDC/hZPAC", "Time Correlation ZP; ZPA Time (ns) ; ZPC Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("All/ZDC/hZEM12", "Time Correlation ZEM; ZEM1 Time (ns) ; ZEM2 Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});

      registry.add("All/ZDC/ZNAamp", "ZNA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZNCamp", "ZNC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZNACCorr", "ZDC amp correlation; ZNA Amplitude; ZNC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/ZDC/ZPAamp", "ZPA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZPCamp", "ZPC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZPACCorr", "ZDC amp correlation; ZPA Amplitude; ZPC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      registry.add("All/ZDC/ZEM1amp", "ZEM1 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZEM2amp", "ZEM2 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("All/ZDC/ZEM12Corr", "ZDC amp correlation; ZEM1 Amplitude; ZEM2 Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/ZDC/ZDCACorr", "ZDC amp correlation; ZDCA Amplitude; ZDCC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // Correlation plots
      registry.add("All/FV0T0ACorr", "Correlation FV0 vs FT0; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FV0T0CCorr", "Correlation FV0 vs FT0; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FT0DDACorr", "Correlation FT0 vs FDD; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FT0DDCCorr", "Correlation FT0 vs FDD; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FT0AFDDC", "Correlation FT0 vs FDD; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FT0CFDDA", "Correlation FT0 vs FDD; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FV0AFDDA", "Correlation FV0 vs FDD; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("All/FV0AFDDC", "Correlation FV0 vs FDD; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
    }

    if (context.mOptions.get<bool>("processHadronic")) {
      registry.add("collHadronic/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
      registry.add("collHadronic/RelBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("collHadronic/trkmult", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("collHadronic/PVTrk", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("collHadronic/etapt", "Eta vs pT; eta of track; pT of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("collHadronic/dEdxTPC", "TPC signal vs signed track pt; Signed track pt [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("collHadronic/dEdxTOF", "TOF signal vs signed track pt; Signed track pt [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
      registry.add("collHadronic/trketa", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("collHadronic/trkpt", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});

      registry.add("collHadronic/trketaZDC", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("collHadronic/trkptZDC", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("collHadronic/trkmultZDC", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("collHadronic/PVTrkZDC", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("collHadronic/etaptZDC", "Eta vs pT; eta of track; pT of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("collHadronic/dEdxTPCZDC", "TPC signal vs signed track pt; Signed track pt [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("collHadronic/dEdxTOFZDC", "TOF signal vs signed track pt; Signed track pt [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // FV0
      registry.add("collHadronic/FV0/hV0A", "Time FV0A; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
      registry.add("collHadronic/FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});

      // FT0
      registry.add("collHadronic/FT0/hT0AC", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("collHadronic/FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("collHadronic/FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("collHadronic/FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // FDD
      registry.add("collHadronic/FDD/hFDDAC", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("collHadronic/FDD/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("collHadronic/FDD/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("collHadronic/FDD/FDDACCorr", "#FDD amp correlation;FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // ZDC
      registry.add("collHadronic/ZDC/ZNAamp", "ZNA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZNCamp", "ZNC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZNACCorr", "ZDC amp correlation; ZNA Amplitude; ZNC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/ZDC/ZPAamp", "ZPA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZPCamp", "ZPC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZPACCorr", "ZDC amp correlation; ZPA Amplitude; ZPC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      registry.add("collHadronic/ZDC/ZEM1amp", "ZEM1 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZEM2amp", "ZEM2 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("collHadronic/ZDC/ZEM12Corr", "ZDC amp correlation; ZEM1 Amplitude; ZEM2 Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // Correlation plots
      registry.add("collHadronic/FV0T0ACorr", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FV0T0CCorr", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FT0DDACorr", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FT0AFDDC", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FT0CFDDA", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FV0AFDDA", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("collHadronic/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
    }

    if (context.mOptions.get<bool>("processInclusiveA")) {

      registry.add("colInclusiveA/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
      registry.add("colInclusiveA/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("colInclusiveA/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveA/PVTracks", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveA/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("colInclusiveA/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("colInclusiveA/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      registry.add("colInclusiveA/trketa", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("colInclusiveA/trkpt", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("colInclusiveA/trketaZDC", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("colInclusiveA/trkptZDC", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("colInclusiveA/trkmultZDC", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveA/PVTracksZDC", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveA/etaptZDC", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("colInclusiveA/dEdxTPCZDC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("colInclusiveA/dEdxTOFZDC", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // FV0
      registry.add("colInclusiveA/FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      // FT0
      registry.add("colInclusiveA/FT0/hT0AC", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("colInclusiveA/FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveA/FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveA/FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      // FDD
      registry.add("colInclusiveA/FDD/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveA/FDD/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveA/FDD/hFDDAC", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0}, {500, -50.0, 50.0}}});
      registry.add("colInclusiveA/FDD/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // ZDC
      registry.add("colInclusiveA/ZDC/ZNAamp", "ZNA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZNCamp", "ZNC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZNACCorr", "ZDC amp correlation; ZNA Amplitude; ZNC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/ZDC/ZPAamp", "ZPA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZPCamp", "ZPC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZPACCorr", "ZDC amp correlation; ZPA Amplitude; ZPC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      registry.add("colInclusiveA/ZDC/ZEM1amp", "ZEM1 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZEM2amp", "ZEM2 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveA/ZDC/ZEM12Corr", "ZDC amp correlation; ZEM1 Amplitude; ZEM2 Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // Correlation plots
      registry.add("colInclusiveA/FV0T0ACorr", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FV0T0CCorr", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FT0DDACorr", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FT0AFDDC", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FT0CFDDA", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FV0AFDDA", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveA/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
    }

    if (context.mOptions.get<bool>("processInclusiveC")) {
      registry.add("colInclusiveC/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
      registry.add("colInclusiveC/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("colInclusiveC/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveC/PVTracks", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveC/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("colInclusiveC/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("colInclusiveC/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
      registry.add("colInclusiveC/trketa", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("colInclusiveC/trkpt", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});

      registry.add("colInclusiveC/trketaZDC", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("colInclusiveC/trkptZDC", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("colInclusiveC/trkmultZDC", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveC/PVTracksZDC", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("colInclusiveC/etaptZDC", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("colInclusiveC/dEdxTPCZDC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("colInclusiveC/dEdxTOFZDC", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // FV0
      registry.add("colInclusiveC/FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      // FT0
      registry.add("colInclusiveC/FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveC/FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveC/FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      // FDD
      registry.add("colInclusiveC/FDD/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveC/FDD/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("colInclusiveC/FDD/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // ZDC
      registry.add("colInclusiveC/ZDC/ZNAamp", "ZNA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZNCamp", "ZNC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZNACCorr", "ZDC amp correlation; ZNA Amplitude; ZNC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/ZDC/ZPAamp", "ZPA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZPCamp", "ZPC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZPACCorr", "ZDC amp correlation; ZPA Amplitude; ZPC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      registry.add("colInclusiveC/ZDC/ZEM1amp", "ZEM1 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZEM2amp", "ZEM2 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("colInclusiveC/ZDC/ZEM12Corr", "ZDC amp correlation; ZEM1 Amplitude; ZEM2 Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // Correlation plots
      registry.add("colInclusiveC/FV0T0ACorr", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FV0T0CCorr", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FT0DDACorr", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FT0AFDDC", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FT0CFDDA", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FV0AFDDA", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("colInclusiveC/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
    }

    if (context.mOptions.get<bool>("processExclusive")) {
      registry.add("exclusive/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
      registry.add("exclusive/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("exclusive/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("exclusive/PVTracks", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("exclusive/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("exclusive/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("exclusive/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
      registry.add("exclusive/trketa", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("exclusive/trkpt", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});

      registry.add("exclusive/trketaZDC", "Eta distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{80, -2., 2.}}});
      registry.add("exclusive/trkptZDC", "Pt distribution of tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{100, 0., 5.}}});
      registry.add("exclusive/trkmultZDC", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("exclusive/PVTracksZDC", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
      registry.add("exclusive/etaptZDC", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("exclusive/dEdxTPCZDC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("exclusive/dEdxTOFZDC", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // FV0
      registry.add("exclusive/FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});

      // FT0
      registry.add("exclusive/FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("exclusive/FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
      registry.add("exclusive/FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      // FDD
      registry.add("exclusive/FDD/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // ZDC
      registry.add("exclusive/ZDC/ZNAamp", "ZNA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZNCamp", "ZNC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZNACCorr", "ZDC amp correlation; ZNA Amplitude; ZNC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/ZDC/ZPAamp", "ZPA Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZPCamp", "ZPC Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZPACCorr", "ZDC amp correlation; ZPA Amplitude; ZPC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      registry.add("exclusive/ZDC/ZEM1amp", "ZEM1 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZEM2amp", "ZEM2 Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
      registry.add("exclusive/ZDC/ZEM12Corr", "ZDC amp correlation; ZEM1 Amplitude; ZEM2 Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

      // Correlation plots
      registry.add("exclusive/FV0T0ACorr", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FV0T0CCorr", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FT0DDACorr", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FT0AFDDC", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FT0CFDDA", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FV0AFDDA", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      registry.add("exclusive/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
    }
  }

  //...............................................................................................................
  void processMain(CC const& collision, BCs const& /*bct0s*/, TCs const& tracks, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/, aod::Zdcs& /*zdcs*/, aod::V0s const& /*v0s*/)
  {
    uint64_t bcnum = 0;
    LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
    float totAmplitudeA = 0;
    float totAmplitudeC = 0;
    float totalAmplitudefv0 = 0;
    float totAmpFddA = 0;
    float totAmpFddC = 0; // auto FITA = 0; auto FITC = 0;
    auto totAmpZNA = 0;
    auto totAmpZNC = 0;
    auto totAmpZPA = 0;
    auto totAmpZPC = 0;
    auto totAmpZEM1 = 0;
    auto totAmpZEM2 = 0; // auto ZDCA = 0; auto ZDCC = 0;

    registry.get<TH1>(HIST("All/Stat"))->Fill(0.);

    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }

    registry.get<TH1>(HIST("All/AllRelBC"))->Fill(bcnum, 1.);
    if (!collision.has_foundBC())
      return;

    registry.get<TH1>(HIST("All/Stat"))->Fill(1.);

    if (collision.has_foundFT0()) {

      auto ft0 = collision.foundFT0();
      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
        registry.get<TH2>(HIST("All/FT0/FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
      }

      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
        registry.get<TH2>(HIST("All/FT0/FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
      }

      for (auto ampa : ft0.amplitudeA()) {
        totAmplitudeA += ampa;
      }

      for (auto ampc : ft0.amplitudeC()) {
        totAmplitudeC += ampc;
      }

      registry.get<TH1>(HIST("All/Stat"))->Fill(2.);
      registry.get<TH2>(HIST("All/FT0/hT0AC"))->Fill(ft0.timeA(), ft0.timeC());

    } else {
      if (!collision.has_foundFT0()) {
        totAmplitudeA = 0;
        totAmplitudeC = 0;
      }
    } // ends FT0

    // FV0 information
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      registry.get<TH1>(HIST("All/Stat"))->Fill(3.);
      registry.get<TH1>(HIST("All/FV0/hV0A"))->Fill(fv0.time());

      for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
        registry.get<TH2>(HIST("All/FV0/FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
      }

      for (auto ampfv0a : fv0.amplitude()) {
        totalAmplitudefv0 += ampfv0a;
      }

    } else {
      if (!collision.has_foundFV0()) {
        totalAmplitudefv0 = 0;
      }
    }

    // FDD information
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("All/Stat"))->Fill(4.);
      registry.get<TH2>(HIST("All/FDD/hFDDAC"))->Fill(fdd.timeA(), fdd.timeC());

      // side A
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("All/FDD/FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
      }

      // side C
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("All/FDD/FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
      }

      for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
      }

      for (auto ampfddc : fdd.chargeC()) {
        totAmpFddC += ampfddc;
      }
    } else {
      if (!collision.has_foundFDD()) {
        totAmpFddA = 0;
        totAmpFddC = 0;
      }
    } // fdd

    // ZDC information
    if (collision.has_foundZDC()) {
      auto zdc = collision.foundZDC();
      registry.get<TH2>(HIST("All/ZDC/hZNAC"))->Fill(zdc.timeZNA(), zdc.timeZNC());
      registry.get<TH2>(HIST("All/ZDC/hZPAC"))->Fill(zdc.timeZPA(), zdc.timeZPC());
      registry.get<TH2>(HIST("All/ZDC/hZEM12"))->Fill(zdc.timeZEM1(), zdc.timeZEM2());

      totAmpZNA = zdc.amplitudeZNA();
      totAmpZNC = zdc.amplitudeZNC();
      totAmpZPA = zdc.amplitudeZPA();
      totAmpZPC = zdc.amplitudeZPC();
      totAmpZEM1 = zdc.amplitudeZEM1();
      totAmpZEM2 = zdc.amplitudeZEM2();
    }

    auto FITA = totalAmplitudefv0 > 0. || totAmplitudeA > 0. || totAmpFddA > 0.;
    auto FITC = totAmplitudeC > 0. || totAmpFddC > 0.;
    auto ZDCA = (totAmpZNA > 0. || totAmpZPA > 0. || totAmpZEM1 > 0. || totAmpZEM2 > 0.);
    auto ZDCC = (totAmpZNC > 0. || totAmpZPC > 0.);

    registry.get<TH1>(HIST("All/RelativeBC"))->Fill(bcnum, 1.);
    registry.get<TH1>(HIST("All/FITAamp"))->Fill(FITA);
    registry.get<TH1>(HIST("All/FITCamp"))->Fill(FITC);
    registry.get<TH2>(HIST("All/FITACCorr"))->Fill(FITA, FITC);

    registry.get<TH1>(HIST("All/ZDC/ZNAamp"))->Fill(totAmpZNA);
    registry.get<TH1>(HIST("All/ZDC/ZNCamp"))->Fill(totAmpZNC);
    registry.get<TH2>(HIST("All/ZDC/ZNACCorr"))->Fill(totAmpZNA, totAmpZNC);
    registry.get<TH2>(HIST("All/ZDC/ZDCACorr"))->Fill(ZDCA, ZDCC);

    registry.get<TH1>(HIST("All/ZDC/ZPAamp"))->Fill(totAmpZPA);
    registry.get<TH1>(HIST("All/ZDC/ZPCamp"))->Fill(totAmpZPC);
    registry.get<TH2>(HIST("All/ZDC/ZPACCorr"))->Fill(totAmpZPA, totAmpZPC);

    registry.get<TH1>(HIST("All/ZDC/ZEM1amp"))->Fill(totAmpZEM1);
    registry.get<TH1>(HIST("All/ZDC/ZEM2amp"))->Fill(totAmpZEM2);
    registry.get<TH2>(HIST("All/ZDC/ZEM12Corr"))->Fill(totAmpZEM1, totAmpZEM2);

    registry.get<TH1>(HIST("All/FV0/FV0Amp"))->Fill(totalAmplitudefv0);
    registry.get<TH1>(HIST("All/FT0/FT0Aamp"))->Fill(totAmplitudeA);
    registry.get<TH1>(HIST("All/FT0/FT0Camp"))->Fill(totAmplitudeC);
    registry.get<TH2>(HIST("All/FT0/FT0ACCorr"))->Fill(totAmplitudeA, totAmplitudeC);
    registry.get<TH1>(HIST("All/FDD/FDDAamp"))->Fill(totAmpFddA);
    registry.get<TH1>(HIST("All/FDD/FDDCamp"))->Fill(totAmpFddC);
    registry.get<TH2>(HIST("All/FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);

    // Correlation FV0 vs FT0
    registry.get<TH2>(HIST("All/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
    registry.get<TH2>(HIST("All/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);

    // Correlation FDD vs FT0
    registry.get<TH2>(HIST("All/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
    registry.get<TH2>(HIST("All/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
    registry.get<TH2>(HIST("All/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
    registry.get<TH2>(HIST("All/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);

    // Correlation FDD vs FV0
    registry.get<TH2>(HIST("All/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
    registry.get<TH2>(HIST("All/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);
    //   }

    // PV contributors
    int nPVcont = 0;
    int nCount = 0;
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.isPVContributor()) {
          nPVcont++;
          registry.get<TH1>(HIST("All/trketa"))->Fill(trk.eta());
          registry.get<TH1>(HIST("All/trkpt"))->Fill(trk.pt());
        }
        nCount++;
        registry.get<TH2>(HIST("All/etapt"))->Fill(trk.eta(), trk.pt(), 1.);
        registry.get<TH2>(HIST("All/dEdxTPC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());

        if (trk.hasTOF()) {
          registry.get<TH2>(HIST("All/dEdxTOF"))->Fill(trk.p() / trk.sign(), trk.beta());
        }
      }
    } // track loop

    registry.get<TH1>(HIST("All/trkmultiplicity"))->Fill(nCount); // all tracks
    registry.get<TH1>(HIST("All/PVTracks"))->Fill(nPVcont);       // PVtracks

    if (collision.has_foundFT0() || collision.has_foundFDD() || collision.has_foundFV0() || collision.has_foundZDC()) {
      registry.get<TH1>(HIST("All/PVCFIT"))->Fill(nPVcont);
    }
  }

  PROCESS_SWITCH(FITtest, processMain, "Process Main", true);
  //...............................................................................................................
  void processHadronic(CC const& collision, BCs const& /*bct0s*/, TCs const& tracks, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/, aod::Zdcs& /*zdcs*/, aod::V0s const& /*v0s*/)
  {
    LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());

    float totAmplitudeA = 0;
    float totAmplitudeC = 0;
    float totalAmplitudefv0 = 0;
    float totAmpFddA = 0;
    float totAmpFddC = 0;
    auto totAmpZNA = 0;
    auto totAmpZNC = 0;
    auto totAmpZPA = 0;
    auto totAmpZPC = 0;
    auto totAmpZEM1 = 0;
    auto totAmpZEM2 = 0;

    uint64_t bcnum = 0;
    registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(0.);

    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }

    if (!collision.has_foundBC())
      return;
    registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(1.);

    if (collision.numContrib() < 10)
      return;
    registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(2.);

    // FIT signal
    if (collision.has_foundFT0()) {
      registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(3.);
      auto ft0 = collision.foundFT0();

      for (auto ampa : ft0.amplitudeA()) {
        totAmplitudeA += ampa;
      }

      for (auto ampc : ft0.amplitudeC()) {
        totAmplitudeC += ampc;
      }
      registry.get<TH2>(HIST("collHadronic/FT0/hT0AC"))->Fill(ft0.timeA(), ft0.timeC());

    } else {
      if (!collision.has_foundFT0()) {
        totAmplitudeA = totAmplitudeC = -999;
      }
    } // ends FT0 collsion

    // FV0 information
    if (collision.has_foundFV0()) {
      registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(4.);

      auto fv0 = collision.foundFV0();
      registry.get<TH1>(HIST("collHadronic/FV0/hV0A"))->Fill(fv0.time());

      for (auto ampfv0a : fv0.amplitude()) {
        totalAmplitudefv0 += ampfv0a;
      }

    } else {
      if (!collision.has_foundFV0()) {
        totalAmplitudefv0 = -999;
      }
    } // FV0 collisions

    // FDD information
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(5.);
      registry.get<TH2>(HIST("collHadronic/FDD/hFDDAC"))->Fill(fdd.timeA(), fdd.timeC());

      for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
      }

      for (auto ampfddc : fdd.chargeC()) {
        totAmpFddC += ampfddc;
      }
    } else {
      if (!collision.has_foundFDD()) {
        totAmpFddA = totAmpFddC = -999;
      }
    } // fdd

    // ZDC information
    if (collision.has_foundZDC()) {
      auto zdc = collision.foundZDC();
      totAmpZNA = zdc.amplitudeZNA();
      totAmpZNC = zdc.amplitudeZNC();
      totAmpZPA = zdc.amplitudeZPA();
      totAmpZPC = zdc.amplitudeZPC();
      totAmpZEM1 = zdc.amplitudeZEM1();
      totAmpZEM2 = zdc.amplitudeZEM2();
    }

    auto FITA = totalAmplitudefv0 > 0. || totAmplitudeA > 0. || totAmpFddA > 0.;
    auto FITC = totAmplitudeC > 0. || totAmpFddC > 0.;
    auto ZDCA = (totAmpZNA > 0. || totAmpZPA > 0. || totAmpZEM1 > 0. || totAmpZEM2 > 0.);
    auto ZDCC = (totAmpZNC > 0. || totAmpZPC > 0.);

    if (!ZDCA || !ZDCC)
      return; // investigate events with ZDC signal on both sides

    registry.get<TH1>(HIST("collHadronic/Stat"))->Fill(6.);
    registry.get<TH1>(HIST("collHadronic/RelBC"))->Fill(bcnum, 1.);

    registry.get<TH1>(HIST("collHadronic/FV0/FV0Amp"))->Fill(totalAmplitudefv0);
    registry.get<TH1>(HIST("collHadronic/FT0/FT0Aamp"))->Fill(totAmplitudeA);
    registry.get<TH1>(HIST("collHadronic/FT0/FT0Camp"))->Fill(totAmplitudeC);
    registry.get<TH2>(HIST("collHadronic/FT0/FT0ACCorr"))->Fill(totAmplitudeA, totAmplitudeC);
    registry.get<TH1>(HIST("collHadronic/FDD/FDDAamp"))->Fill(totAmpFddA);
    registry.get<TH1>(HIST("collHadronic/FDD/FDDCamp"))->Fill(totAmpFddC);
    registry.get<TH2>(HIST("collHadronic/FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);

    // Correlation FV0 vs FT0
    registry.get<TH2>(HIST("collHadronic/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
    registry.get<TH2>(HIST("collHadronic/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);

    // Correlation FDD vs FT0
    registry.get<TH2>(HIST("collHadronic/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
    registry.get<TH2>(HIST("collHadronic/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
    registry.get<TH2>(HIST("collHadronic/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
    registry.get<TH2>(HIST("collHadronic/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);

    // Correlation FDD vs FV0
    registry.get<TH2>(HIST("collHadronic/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
    registry.get<TH2>(HIST("collHadronic/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);

    // PV contributors
    int nPVcont = 0;
    int nCont = 0;
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() > 1) {
          if (trk.isPVContributor()) {
            nPVcont++;
            registry.get<TH1>(HIST("collHadronic/trketa"))->Fill(trk.eta());
            registry.get<TH1>(HIST("collHadronic/trkpt"))->Fill(trk.pt());
          }
          nCont++;
          registry.get<TH2>(HIST("collHadronic/etapt"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("collHadronic/dEdxTPC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("collHadronic/dEdxTOF"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // tracks

    registry.get<TH1>(HIST("collHadronic/PVTrk"))->Fill(nPVcont);
    registry.get<TH1>(HIST("collHadronic/trkmult"))->Fill(nCont);

    if (!FITA || !FITC)
      return; // investigate events with FIT signal on both sides

    registry.get<TH1>(HIST("collHadronic/ZDC/ZNAamp"))->Fill(totAmpZNA);
    registry.get<TH1>(HIST("collHadronic/ZDC/ZNCamp"))->Fill(totAmpZNC);
    registry.get<TH2>(HIST("collHadronic/ZDC/ZNACCorr"))->Fill(totAmpZNA, totAmpZNC);

    registry.get<TH1>(HIST("collHadronic/ZDC/ZPAamp"))->Fill(totAmpZPA);
    registry.get<TH1>(HIST("collHadronic/ZDC/ZPCamp"))->Fill(totAmpZPC);
    registry.get<TH2>(HIST("collHadronic/ZDC/ZPACCorr"))->Fill(totAmpZPA, totAmpZPC);

    registry.get<TH1>(HIST("collHadronic/ZDC/ZEM1amp"))->Fill(totAmpZEM1);
    registry.get<TH1>(HIST("collHadronic/ZDC/ZEM2amp"))->Fill(totAmpZEM2);
    registry.get<TH2>(HIST("collHadronic/ZDC/ZEM12Corr"))->Fill(totAmpZEM1, totAmpZEM2);

    // PV contributors
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() > 1) {
          if (trk.isPVContributor()) {
            registry.get<TH1>(HIST("collHadronic/trketaZDC"))->Fill(trk.eta());
            registry.get<TH1>(HIST("collHadronic/trkptZDC"))->Fill(trk.pt());
          }
          registry.get<TH2>(HIST("collHadronic/etaptZDC"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("collHadronic/dEdxTPCZDC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("collHadronic/dEdxTOFZDC"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // tracks

    registry.get<TH1>(HIST("collHadronic/PVTrkZDC"))->Fill(nPVcont);
    registry.get<TH1>(HIST("collHadronic/trkmultZDC"))->Fill(nCont);
  }

  PROCESS_SWITCH(FITtest, processHadronic, "Process for hadroniclike events", true);
  //...............................................................................................................
  void processInclusiveA(CC const& collision, BCs const& /*bct0s*/, TCs const& tracks, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/, aod::Zdcs& /*zdcs*/, aod::V0s const& /*v0s*/)
  {
    uint64_t bcnum = 0;
    float totAmplitudeA = 0;
    float totAmplitudeC = 0;
    float totalAmplitudefv0 = 0;
    float totAmpFddA = 0;
    float totAmpFddC = 0;
    auto totAmpZNA = 0;
    auto totAmpZNC = 0;
    auto totAmpZPA = 0;
    auto totAmpZPC = 0;
    auto totAmpZEM1 = 0;
    auto totAmpZEM2 = 0;

    LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
    registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(0.);

    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }

    if (!collision.has_foundBC())
      return;
    registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(1.);

    if (collision.numContrib() > 50)
      return;
    registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(2.);

    // FT0 information
    if (collision.has_foundFT0()) {
      registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(3.);
      auto ft0 = collision.foundFT0();

      for (auto ampa : ft0.amplitudeA()) {
        totAmplitudeA += ampa;
      }

      for (auto ampc : ft0.amplitudeC()) {
        totAmplitudeC += ampc;
      }
      registry.get<TH2>(HIST("colInclusiveA/FT0/hT0AC"))->Fill(ft0.timeA(), ft0.timeC());

    } else {
      if (!collision.has_foundFT0()) {
        totAmplitudeA = 0;
        totAmplitudeC = 0;
      }
    } // ends collsion

    // FV0 information
    if (collision.has_foundFV0()) {
      registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(4.);

      auto fv0 = collision.foundFV0();

      for (auto ampfv0a : fv0.amplitude()) {
        totalAmplitudefv0 += ampfv0a;
      }

    } else {
      if (!collision.has_foundFV0()) {
        totalAmplitudefv0 = 0;
      }
    }

    // FDD information
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(5.);
      registry.get<TH2>(HIST("colInclusiveA/FDD/hFDDAC"))->Fill(fdd.timeA(), fdd.timeC());

      for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
      }

      for (auto ampfddc : fdd.chargeC()) {
        totAmpFddC += ampfddc;
      }
    } else {
      if (!collision.has_foundFDD()) {
        totAmpFddA = 0;
        totAmpFddC = 0;
      }
    } // fdd

    // ZDC
    if (collision.has_foundZDC()) {
      auto zdc = collision.foundZDC();
      totAmpZNA = zdc.amplitudeZNA();
      totAmpZNC = zdc.amplitudeZNC();
      totAmpZPA = zdc.amplitudeZPA();
      totAmpZPC = zdc.amplitudeZPC();
      totAmpZEM1 = zdc.amplitudeZEM1();
      totAmpZEM2 = zdc.amplitudeZEM2();
    }

    auto FITA = (totalAmplitudefv0 > 0. || totAmplitudeA > 0. || totAmpFddA > 0.);
    auto FITC = (totAmplitudeC > 0. || totAmpFddC > 0.);
    auto ZDCA = (totAmpZNA > 0. || totAmpZPA > 0. || totAmpZEM1 > 0. || totAmpZEM2 > 0.);
    auto ZDCC = (totAmpZNC > 0. || totAmpZPC > 0.);

    if (!ZDCA || ZDCC)
      return;
    registry.get<TH1>(HIST("colInclusiveA/RelativeBC"))->Fill(bcnum, 1.);
    registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(6.);

    registry.get<TH1>(HIST("colInclusiveA/FV0/FV0Amp"))->Fill(totalAmplitudefv0);
    registry.get<TH1>(HIST("colInclusiveA/FT0/FT0Aamp"))->Fill(totAmplitudeA);
    registry.get<TH1>(HIST("colInclusiveA/FT0/FT0Camp"))->Fill(totAmplitudeC);
    registry.get<TH2>(HIST("colInclusiveA/FT0/FT0ACCorr"))->Fill(totAmplitudeA, totAmplitudeC);
    registry.get<TH2>(HIST("colInclusiveA/FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);
    registry.get<TH1>(HIST("colInclusiveA/FDD/FDDAamp"))->Fill(totAmpFddA);
    registry.get<TH1>(HIST("colInclusiveA/FDD/FDDCamp"))->Fill(totAmpFddC);

    // Correlation FV0 vs FT0
    registry.get<TH2>(HIST("colInclusiveA/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
    registry.get<TH2>(HIST("colInclusiveA/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);

    // Correlation FDD vs FT0
    registry.get<TH2>(HIST("colInclusiveA/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveA/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
    registry.get<TH2>(HIST("colInclusiveA/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveA/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);

    // Correlation FDD vs FV0
    registry.get<TH2>(HIST("colInclusiveA/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveA/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);

    // PV contributors
    int nPVcont = 0;
    int ntrks = 0;
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            nPVcont++;
            registry.get<TH1>(HIST("colInclusiveA/trketa"))->Fill(trk.eta());
            registry.get<TH1>(HIST("colInclusiveA/trkpt"))->Fill(trk.pt());
          }

          ntrks++;
          registry.get<TH2>(HIST("colInclusiveA/etapt"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("colInclusiveA/dEdxTPC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("colInclusiveA/dEdxTOF"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // trk loop
    registry.get<TH1>(HIST("colInclusiveA/PVTracks"))->Fill(nPVcont);
    registry.get<TH1>(HIST("colInclusiveA/trkmultiplicity"))->Fill(ntrks);

    if (!FITA || FITC)
      return;
    registry.get<TH1>(HIST("colInclusiveA/Stat"))->Fill(7.);

    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZNAamp"))->Fill(totAmpZNA);
    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZNCamp"))->Fill(totAmpZNC);
    registry.get<TH2>(HIST("colInclusiveA/ZDC/ZNACCorr"))->Fill(totAmpZNA, totAmpZNC);

    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZPAamp"))->Fill(totAmpZPA);
    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZPCamp"))->Fill(totAmpZPC);
    registry.get<TH2>(HIST("colInclusiveA/ZDC/ZPACCorr"))->Fill(totAmpZPA, totAmpZPC);

    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZEM1amp"))->Fill(totAmpZEM1);
    registry.get<TH1>(HIST("colInclusiveA/ZDC/ZEM2amp"))->Fill(totAmpZEM2);
    registry.get<TH2>(HIST("colInclusiveA/ZDC/ZEM12Corr"))->Fill(totAmpZEM1, totAmpZEM2);

    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            registry.get<TH1>(HIST("colInclusiveA/trketaZDC"))->Fill(trk.eta());
            registry.get<TH1>(HIST("colInclusiveA/trkptZDC"))->Fill(trk.pt());
          }
          registry.get<TH2>(HIST("colInclusiveA/etaptZDC"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("colInclusiveA/dEdxTPCZDC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("colInclusiveA/dEdxTOFZDC"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // trk loop

    registry.get<TH1>(HIST("colInclusiveA/PVTracksZDC"))->Fill(nPVcont);
    registry.get<TH1>(HIST("colInclusiveA/trkmultZDC"))->Fill(ntrks);
  }

  PROCESS_SWITCH(FITtest, processInclusiveA, "Process Inclusive veto A side", true);
  //..................................................................................................................................
  void processInclusiveC(CC const& collision, BCs const& /*bct0s*/, TCs const& tracks, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/, aod::Zdcs& /*zdcs*/, aod::V0s const& /*v0s*/)
  {
    uint64_t bcnum = 0;
    float totAmplitudeA = 0;
    float totAmplitudeC = 0;
    float totalAmplitudefv0 = 0;
    float totAmpFddA = 0;
    float totAmpFddC = 0;
    auto totAmpZNA = 0;
    auto totAmpZNC = 0;
    auto totAmpZPA = 0;
    auto totAmpZPC = 0;
    auto totAmpZEM1 = 0;
    auto totAmpZEM2 = 0;

    LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
    registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(0.);

    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }

    if (!collision.has_foundBC())
      return;
    registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(1.);

    if (collision.numContrib() > 50)
      return;
    registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(2.);

    // FT0 information
    if (collision.has_foundFT0()) {
      registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(3.);
      auto ft0 = collision.foundFT0();

      for (auto ampa : ft0.amplitudeA()) {
        totAmplitudeA += ampa;
      }

      for (auto ampc : ft0.amplitudeC()) {
        totAmplitudeC += ampc;
      }

    } else {
      if (!collision.has_foundFT0()) {
        totAmplitudeA = 0;
        totAmplitudeC = 0;
      }
    } // ends collsion

    // FV0 information
    if (collision.has_foundFV0()) {
      registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(4.);

      auto fv0 = collision.foundFV0();
      for (auto ampfv0a : fv0.amplitude()) {
        totalAmplitudefv0 += ampfv0a;
      }

    } else {
      if (!collision.has_foundFV0()) {
        totalAmplitudefv0 = 0;
      }
    }

    // FDD information
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(5.);

      for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
      }

      for (auto ampfddc : fdd.chargeC()) {
        totAmpFddC += ampfddc;
      }
    } else {
      if (!collision.has_foundFDD()) {
        totAmpFddA = 0;
        totAmpFddC = 0;
      }
    } // fdd

    // ZDC
    if (collision.has_foundZDC()) {
      auto zdc = collision.foundZDC();
      totAmpZNA = zdc.amplitudeZNA();
      totAmpZNC = zdc.amplitudeZNC();
      totAmpZPA = zdc.amplitudeZPA();
      totAmpZPC = zdc.amplitudeZPC();
      totAmpZEM1 = zdc.amplitudeZEM1();
      totAmpZEM2 = zdc.amplitudeZEM2();
    }

    auto FITA = totalAmplitudefv0 > 0. || totAmplitudeA > 0. || totAmpFddA > 0.;
    auto FITC = totAmplitudeC > 0. || totAmpFddC > 0.;
    auto ZDCA = (totAmpZNA > 0. || totAmpZPA > 0. || totAmpZEM1 > 0. || totAmpZEM2 > 0.);
    auto ZDCC = (totAmpZNC > 0. || totAmpZPC > 0.);

    if (ZDCA || !ZDCC)
      return;
    registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(6.);
    registry.get<TH1>(HIST("colInclusiveC/RelativeBC"))->Fill(bcnum, 1.);

    registry.get<TH1>(HIST("colInclusiveC/FV0/FV0Amp"))->Fill(totalAmplitudefv0);
    registry.get<TH1>(HIST("colInclusiveC/FT0/FT0Aamp"))->Fill(totAmplitudeA);
    registry.get<TH1>(HIST("colInclusiveC/FT0/FT0Camp"))->Fill(totAmplitudeC);
    registry.get<TH2>(HIST("colInclusiveC/FT0/FT0ACCorr"))->Fill(totAmplitudeA, totAmplitudeC);
    registry.get<TH2>(HIST("colInclusiveC/FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);
    registry.get<TH1>(HIST("colInclusiveC/FDD/FDDAamp"))->Fill(totAmpFddA);
    registry.get<TH1>(HIST("colInclusiveC/FDD/FDDCamp"))->Fill(totAmpFddC);

    // Correlation FV0 vs FT0
    registry.get<TH2>(HIST("colInclusiveC/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
    registry.get<TH2>(HIST("colInclusiveC/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);

    // Correlation FDD vs FT0
    registry.get<TH2>(HIST("colInclusiveC/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveC/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
    registry.get<TH2>(HIST("colInclusiveC/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveC/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);

    // Correlation FDD vs FV0
    registry.get<TH2>(HIST("colInclusiveC/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
    registry.get<TH2>(HIST("colInclusiveC/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);

    // PV contributors
    int nPVcont = 0;
    int ntrks = 0;
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            nPVcont++;
            registry.get<TH1>(HIST("colInclusiveC/trketa"))->Fill(trk.eta());
            registry.get<TH1>(HIST("colInclusiveC/trkpt"))->Fill(trk.pt());
          }
          ntrks++;
          registry.get<TH2>(HIST("colInclusiveC/etapt"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("colInclusiveC/dEdxTPC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("colInclusiveC/dEdxTOF"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    }
    registry.get<TH1>(HIST("colInclusiveC/PVTracks"))->Fill(nPVcont);
    registry.get<TH1>(HIST("colInclusiveC/trkmultiplicity"))->Fill(ntrks);

    if (FITA || !FITC)
      return;
    registry.get<TH1>(HIST("colInclusiveC/Stat"))->Fill(7.);

    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZNAamp"))->Fill(totAmpZNA);
    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZNCamp"))->Fill(totAmpZNC);
    registry.get<TH2>(HIST("colInclusiveC/ZDC/ZNACCorr"))->Fill(totAmpZNA, totAmpZNC);

    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZPAamp"))->Fill(totAmpZPA);
    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZPCamp"))->Fill(totAmpZPC);
    registry.get<TH2>(HIST("colInclusiveC/ZDC/ZPACCorr"))->Fill(totAmpZPA, totAmpZPC);

    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZEM1amp"))->Fill(totAmpZEM1);
    registry.get<TH1>(HIST("colInclusiveC/ZDC/ZEM2amp"))->Fill(totAmpZEM2);
    registry.get<TH2>(HIST("colInclusiveC/ZDC/ZEM12Corr"))->Fill(totAmpZEM1, totAmpZEM2);

    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            registry.get<TH1>(HIST("colInclusiveC/trketaZDC"))->Fill(trk.eta());
            registry.get<TH1>(HIST("colInclusiveC/trkptZDC"))->Fill(trk.pt());
          }
          registry.get<TH2>(HIST("colInclusiveC/etaptZDC"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("colInclusiveC/dEdxTPCZDC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("colInclusiveC/dEdxTOFZDC"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // trk loop

    registry.get<TH1>(HIST("colInclusiveC/PVTracksZDC"))->Fill(nPVcont);
    registry.get<TH1>(HIST("colInclusiveC/trkmultZDC"))->Fill(ntrks);
  }

  PROCESS_SWITCH(FITtest, processInclusiveC, "Process Inclusive veto C side", true);
  //..................................................................................................................................
  void processExclusive(CC const& collision, BCs const& /*bct0s*/, TCs const& tracks, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/, aod::Zdcs& /*zdcs*/, aod::V0s const& /*v0s*/)
  {
    float totAmplitudeA = 0;
    float totAmplitudeC = 0;
    float totalAmplitudefv0 = 0;
    float totAmpFddA = 0;
    float totAmpFddC = 0;
    auto totAmpZNA = 0;
    auto totAmpZNC = 0;
    auto totAmpZPA = 0;
    auto totAmpZPC = 0;
    auto totAmpZEM1 = 0;
    auto totAmpZEM2 = 0;

    uint64_t bcnum = 0;
    // auto bc = collision.foundBC_as<BCs>();
    registry.get<TH1>(HIST("exclusive/Stat"))->Fill(0.);
    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }

    if (!collision.has_foundBC())
      return;
    if (collision.numContrib() > 10)
      return;

    registry.get<TH1>(HIST("exclusive/Stat"))->Fill(1.);
    registry.get<TH1>(HIST("exclusive/trkmultiplicity"))->Fill(tracks.size(), 1.);

    if (collision.has_foundFT0()) {

      registry.get<TH1>(HIST("exclusive/Stat"))->Fill(2.);
      auto ft0 = collision.foundFT0();

      for (auto ampa : ft0.amplitudeA()) {
        totAmplitudeA += ampa;
      }

      for (auto ampc : ft0.amplitudeC()) {
        totAmplitudeC += ampc;
      }

    } else {
      if (!collision.has_foundFT0()) {
        totAmplitudeA = 0;
        totAmplitudeC = 0;
      }
    } // ends collsion

    // FV0 information
    if (collision.has_foundFV0()) {
      registry.get<TH1>(HIST("exclusive/Stat"))->Fill(3.);

      auto fv0 = collision.foundFV0();
      for (auto ampfv0a : fv0.amplitude()) {
        totalAmplitudefv0 += ampfv0a;
      }

    } else {
      if (!collision.has_foundFV0()) {
        totalAmplitudefv0 = 0;
      }
    }

    // FDD information
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("exclusive/Stat"))->Fill(4.);

      for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
      }

      for (auto ampfddc : fdd.chargeC()) {
        totAmpFddC += ampfddc;
      }
    } else {
      if (!collision.has_foundFDD()) {
        totAmpFddA = 0;
        totAmpFddC = 0;
      }
    } // fdd

    // ZDC
    if (collision.has_foundZDC()) {
      auto zdc = collision.foundZDC();
      totAmpZNA = zdc.amplitudeZNA();
      totAmpZNC = zdc.amplitudeZNC();
      totAmpZPA = zdc.amplitudeZPA();
      totAmpZPC = zdc.amplitudeZPC();
      totAmpZEM1 = zdc.amplitudeZEM1();
      totAmpZEM2 = zdc.amplitudeZEM2();
    }

    auto FITA = totalAmplitudefv0 > 0. || totAmplitudeA > 0. || totAmpFddA > 0.;
    auto FITC = totAmplitudeC > 0. || totAmpFddC > 0.;
    auto ZDCA = (totAmpZNA > 0. || totAmpZPA > 0. || totAmpZEM1 > 0. || totAmpZEM2 > 0.);
    auto ZDCC = (totAmpZNC > 0. || totAmpZPC > 0.);

    if (ZDCA || ZDCC)
      return;

    registry.get<TH1>(HIST("exclusive/FT0/FT0Aamp"))->Fill(totAmplitudeA);
    registry.get<TH1>(HIST("exclusive/FT0/FT0Camp"))->Fill(totAmplitudeC);
    registry.get<TH2>(HIST("exclusive/FT0/FT0ACCorr"))->Fill(totAmplitudeA, totAmplitudeC);

    registry.get<TH1>(HIST("exclusive/FV0/FV0Amp"))->Fill(totalAmplitudefv0);
    registry.get<TH2>(HIST("exclusive/FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);

    registry.get<TH2>(HIST("exclusive/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
    registry.get<TH2>(HIST("exclusive/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);

    registry.get<TH2>(HIST("exclusive/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
    registry.get<TH2>(HIST("exclusive/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
    registry.get<TH2>(HIST("exclusive/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
    registry.get<TH2>(HIST("exclusive/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);
    registry.get<TH2>(HIST("exclusive/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
    registry.get<TH2>(HIST("exclusive/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);

    // PV contributors
    int nPVcont = 0;
    int ntrks = 0;
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            nPVcont++;
            registry.get<TH1>(HIST("exclusive/trketa"))->Fill(trk.eta());
            registry.get<TH1>(HIST("exclusive/trkpt"))->Fill(trk.pt());
          }
          ntrks++;
          registry.get<TH2>(HIST("exclusive/etapt"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("exclusive/dEdxTPC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("exclusive/dEdxTOF"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    }

    registry.get<TH1>(HIST("exclusive/PVTracks"))->Fill(nPVcont);

    if (FITA || FITC)
      return;
    registry.get<TH1>(HIST("exclusive/ZDC/ZNAamp"))->Fill(totAmpZNA);
    registry.get<TH1>(HIST("exclusive/ZDC/ZNCamp"))->Fill(totAmpZNC);
    registry.get<TH2>(HIST("exclusive/ZDC/ZNACCorr"))->Fill(totAmpZNA, totAmpZNC);

    registry.get<TH1>(HIST("exclusive/ZDC/ZPAamp"))->Fill(totAmpZPA);
    registry.get<TH1>(HIST("exclusive/ZDC/ZPCamp"))->Fill(totAmpZPC);
    registry.get<TH2>(HIST("exclusive/ZDC/ZPACCorr"))->Fill(totAmpZPA, totAmpZPC);

    registry.get<TH1>(HIST("exclusive/ZDC/ZEM1amp"))->Fill(totAmpZEM1);
    registry.get<TH1>(HIST("exclusive/ZDC/ZEM2amp"))->Fill(totAmpZEM2);
    registry.get<TH2>(HIST("exclusive/ZDC/ZEM12Corr"))->Fill(totAmpZEM1, totAmpZEM2);

    // PV contributors
    for (auto const& trk : tracks) {
      if (trk.eta() > -1.5 && trk.eta() < 1.5) {
        if (trk.pt() < 10) {
          if (trk.isPVContributor()) {
            registry.get<TH1>(HIST("exclusive/trketaZDC"))->Fill(trk.eta());
            registry.get<TH1>(HIST("exclusive/trkptZDC"))->Fill(trk.pt());
          }
          registry.get<TH2>(HIST("exclusive/etaptZDC"))->Fill(trk.eta(), trk.pt(), 1.);
          registry.get<TH2>(HIST("exclusive/dEdxTPCZDC"))->Fill(trk.tpcInnerParam() / trk.sign(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.get<TH2>(HIST("exclusive/dEdxTOFZDC"))->Fill(trk.p() / trk.sign(), trk.beta());
          }
        }
      }
    } // trk loop

    registry.get<TH1>(HIST("exclusive/PVTracksZDC"))->Fill(nPVcont);
    registry.get<TH1>(HIST("exclusive/RelativeBC"))->Fill(bcnum, 1.);
    registry.get<TH1>(HIST("exclusive/trkmultZDC"))->Fill(ntrks);
  }

  PROCESS_SWITCH(FITtest, processExclusive, "Process exclusiveUPC veto A and C sides", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FITtest>(cfgc, TaskName{"fittest"}),
  };
}
