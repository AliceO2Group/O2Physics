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

/// \file tableMakerMuonMchTrkEfficiency.cxx
/// \brief task to prepare the tables and some plots required for evaluating the muon tracking efficiency uncertainty
///
/// @param muon MCH tracking efficiency table
/// Struct for filling the histos and writing the table needed to compute the
///    muon tracking efficiency in the MCH detector
///
/// \author Zaida Conesa del Valle <zaida.conesa.del.valle@cern.ch>
///

#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include <TLorentzVector.h>
#include "TDatabasePDG.h"
//
#include "Common/DataModel/TrackSelectionTables.h"
//
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"
//
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/DataModel/MchTrkEffTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// default histogram output binning
namespace muon_trk_eff_bins
{
static constexpr int nBinsPt = 24;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
  2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,
  15.0, 18.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// row labels
static const std::vector<std::string> labelsPtTrack{
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11",
  "pT bin 12",
  "pT bin 13",
  "pT bin 14",
  "pT bin 15",
  "pT bin 16",
  "pT bin 17",
  "pT bin 18",
  "pT bin 19",
  "pT bin 20",
  "pT bin 21",
  "pT bin 22",
  "pT bin 23",
  "pT bin 24"};

} // namespace muon_trk_eff_bins
// Constants
static const double muonMass = 0.105658; // in GeV from PDG
static const int muonPDG = 13;

///  muon tracking efficiency  task
struct tableMakerMuonMchTrkEfficiency {

  // Declare tables
  Produces<aod::MchTrkEffBase> rowCandidateBase;
  Produces<aod::MchTrkEffGen> rowCandidateGen;

  /// Configure the task variables
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};                   /// Event selection list
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"}; /// List of muon selections
  ///
  Configurable<int> muonSelType{"muonSelType", -1, "Selected muon track type"}; /// Muon track type to be selected if value >=0 (no selection by default)
  Configurable<double> ptMuonMin{"ptMin", 0., "Lower bound of pT"};             /// Muon minimum pt to be studied
  Configurable<double> etaMuonMin{"etaMin", 2.5, "Lower bound of |eta|"};       /// Muon minimum |eta| to be studied
  Configurable<double> etaMuonMax{"etaMax", 4.0, "Upper bound of |eta|"};       /// Muon maximum |eta| to be studied
  ///
  Configurable<std::vector<double>> binsMuonPt{"binsPt", std::vector<double>{muon_trk_eff_bins::vecBinsPt}, "pT bin limits"}; /// Pt intervals for the histograms
  Configurable<int> nEtaBins{"nEtaBins", 12, "Number of Eta bins"};                                                           /// Number of eta bins for output histograms
  Configurable<int> nPhiBins{"nPhiBins", 6, "Number of Phi bins"};                                                            /// Number of phi bins for output histograms
  Configurable<bool> fillBitMapCorr{"fillCorr", false, "Fill bit map correlation sparse"};                                    /// Boolean to fill or not the THnSparse of correlations

  AnalysisCompositeCut* fEventCut;             //! Event selection cut
  std::vector<AnalysisCompositeCut> fMuonCuts; //! Muon track cuts

  /// Declarations of various short names
  using myEvents = soa::Join<aod::Collisions, aod::EvSels>;
  using myEventsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using myReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
  using myReducedEventsMC = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;

  using myMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
  using myMuonsMC = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
  using myReducedMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
  using myReducedMuonsMC = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;

  // bit maps used for the Fill functions of the VarManager
  constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
  constexpr static uint32_t gkReducedEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
  constexpr static uint32_t gkEventMCFillMap = VarManager::ObjTypes::CollisionMC;
  constexpr static uint32_t gkReducedEventMCFillMap = VarManager::ObjTypes::ReducedEventMC;

  constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
  constexpr static uint32_t gkReducedMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;
  constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

  /// Histogram registry: an object to hold your histograms
  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  ///  Initialize: configure, create specifics
  void init(o2::framework::InitContext&)
  {
    LOGF(debug, "Initialization");

    /// set event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    LOGF(debug, ">> Event cut: %s", eventCutStr.Data());

    /// set muon cuts
    TString cutNamesStr = fConfigMuonCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        LOGF(debug, ">> Muon cut added: %s", objArray->At(icut)->GetName());
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    VarManager::SetDefaultVarNames();

    // define histograms to be added
    LOGF(debug, " Creating histograms");
    const AxisSpec axisEvt{10, -0.5, 9.5, ""};
    const char* elabels[6] = {"all", "selected", "sel >0 muon", "sel >1 muon", "sel >2 muon", "sel >3 muon"};
    registry.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {axisEvt}});
    auto hEvent = registry.get<TH1>(HIST("hEventCounter"));
    for (int i = 0; i < 6; i++)
      hEvent->GetXaxis()->SetBinLabel(i + 1, elabels[i]);

    // define axes to be used
    // and the correspondent histograms
    auto vbins = (std::vector<double>)binsMuonPt;
    const AxisSpec axisEta{190, 2.3, 4.2, "|#eta|"};
    const AxisSpec axisEtaRed{nEtaBins, etaMuonMin, etaMuonMax, "|#eta|"};
    const AxisSpec axisEtaGenRed{nEtaBins, etaMuonMin, etaMuonMax, "|#eta| Gen"};
    const AxisSpec axisPt{vbins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtGen{vbins, "#it{p}_{T} (GeV/#it{c}) Gen"};
    const AxisSpec axisPhi{120, -3.14, 3.14, "#varphi"};
    const AxisSpec axisPhiRed{nPhiBins, -3.14, 3.14, "#varphi"};
    const AxisSpec axisPhiGenRed{nPhiBins, -3.14, 3.14, "#varphi Gen"};
    const AxisSpec axisBitMap{1031, -0.5, 1030.5, "mchbitMap"};
    const int nChHitBins = 4;
    // Number of hits in chamber i = (Nij, Ni0, N0j) correspond to hit on i-j, hit on i not on j, hit on j not on i
    const AxisSpec axisNhits{nChHitBins, 0.5, nChHitBins + 0.5, "isHitInChamber"};
    const AxisSpec axisX{60, -3., 3., "X"};
    const AxisSpec axisY{60, -3., 3., "Y"};

    HistogramConfigSpec defaultNhitsPerChamber({HistType::kTH3F, {{axisNhits}, {axisX}, {axisY}}});
    HistogramConfigSpec defaultMchBitmap({HistType::kTHnF, {axisBitMap, axisEtaRed, axisPt, axisPhiRed}});

    registry.add("hEta", "hEta", {HistType::kTH1F, {axisEta}});
    registry.add("hPt", "hPt", {HistType::kTH1F, {axisPt}});
    registry.add("hPhi", "hPhi", {HistType::kTH1F, {axisPhi}});
    registry.add("hMchBitMap", "hMchBitMap", {HistType::kTH1F, {axisBitMap}});

    // keep also histos per muon selection
    registry.add("selected/hEta", Form("hEta_%s", cutNamesStr.Data()), {HistType::kTH1F, {axisEta}});
    registry.add("selected/hPt", Form("hPt_%s", cutNamesStr.Data()), {HistType::kTH1F, {axisPt}});
    registry.add("selected/hPhi", Form("hPhi_%s", cutNamesStr.Data()), {HistType::kTH1F, {axisPhi}});
    registry.add("selected/hMchBitMap", Form("hMchBitMap_%s", cutNamesStr.Data()), {HistType::kTH1F, {axisBitMap}});

    registry.add("hPtRecPtGen", "hPtRecPtGen", {HistType::kTH2F, {{axisPt}, {axisPtGen}}});
    registry.add("hEtaRecEtaGen", "hEtaRecEtaGen", {HistType::kTH2F, {{axisEtaRed}, {axisEtaGenRed}}});
    registry.add("hPhiRecPhiGen", "hPhiRecPhiGen", {HistType::kTH2F, {{axisPhiRed}, {axisPhiGenRed}}});
    registry.add("selected/hPtRecPtGen", "hPtRecPtGen", {HistType::kTH2F, {{axisPt}, {axisPtGen}}});
    registry.add("selected/hEtaRecEtaGen", "hEtaRecEtaGen", {HistType::kTH2F, {{axisEtaRed}, {axisEtaGenRed}}});
    registry.add("selected/hPhiRecPhiGen", "hPhiRecPhiGen", {HistType::kTH2F, {{axisPhiRed}, {axisPhiGenRed}}});

    registry.add("hEtaPtPhi", "hEtaPtPhi", {HistType::kTH3F, {{axisEtaRed}, {axisPt}, {axisPhiRed}}});
    registry.add("selected/hEtaPtPhi", Form("hEtaPtPhi_%s", cutNamesStr.Data()), {HistType::kTH3F, {{axisEtaRed}, {axisPt}, {axisPhiRed}}});

    for (int i = 0; i < 10; i++) {
      registry.add(Form("selected/hNhitsPerChamber_%i", i), Form("hNhitsPerChamber_%i_%s", i, cutNamesStr.Data()), defaultNhitsPerChamber, false);
    }

    registry.add("selected/hMchBitMapEtaPtPhi", Form("hMchBitMapEtaPtPhi_%s", cutNamesStr.Data()), defaultMchBitmap);

    LOGF(debug, "End of initialization");

  } //! end of Initialize: configure, create specifics

  /// check whether a given chamber has hits
  bool ischamberhit(uint16_t map, int ich)
  { // i = 0..9
    LOGF(debug, " map %i --> %i", map, (map >> ich) & 1);
    return (map >> ich) & 1;
  }

  /// extrapolate tracks to a given r value (spherical coordinates)
  ///   to mimic the (x,y) position in a given chamber
  void extrapolate(TLorentzVector vec, int ich, double& x, double& y)
  { // i = 0..9
    double zposCh[10] = {5, 5, 7, 7, 10, 10, 12.5, 12.5, 14.5, 14.5};
    double theta = vec.Theta();
    double phi = vec.Phi();
    double r = zposCh[ich] / TMath::Cos(theta);
    double myx = r * TMath::Sin(theta) * TMath::Cos(phi);
    double myy = r * TMath::Sin(theta) * TMath::Sin(phi);
    x = myx;
    y = myy;
  }

  /// process to fill histograms
  void FillHistos(double mEta, double mPhi, double mPt, uint16_t mchBitmap, bool isSel)
  {

    /// fill histograms
    registry.fill(HIST("hEta"), mEta);
    registry.fill(HIST("hPt"), mPt);
    registry.fill(HIST("hPhi"), mPhi);
    registry.fill(HIST("hMchBitMap"), mchBitmap);
    registry.fill(HIST("hEtaPtPhi"), mEta, mPt, mPhi);

    /// fill histograms only for selected candidates
    if (isSel) {

      registry.fill(HIST("selected/hEta"), mEta);
      registry.fill(HIST("selected/hPt"), mPt);
      registry.fill(HIST("selected/hPhi"), mPhi);
      registry.fill(HIST("selected/hMchBitMap"), mchBitmap);
      registry.fill(HIST("selected/hEtaPtPhi"), mEta, mPt, mPhi);
      if (fillBitMapCorr)
        registry.fill(HIST("selected/hMchBitMapEtaPtPhi"), mchBitmap, mEta, mPt, mPhi);

      /// Study the Nhit distribution vs X-Y
      const int nChambers = 10;
      bool isNChamberHit[nChambers] = {0};
      double xCh[nChambers], yCh[nChambers];
      TLorentzVector mvector;
      mvector.SetPtEtaPhiM(mPt, mEta, mPhi, muonMass);
      for (int i = 0; i < nChambers; i++) {
        xCh[i] = 0.;
        yCh[i] = 0.;
        isNChamberHit[i] = false;
        isNChamberHit[i] = ischamberhit(mchBitmap, i);
        extrapolate(mvector, i, xCh[i], yCh[i]);
      }

      // Fill histos station 1, Chambers 1 & 2
      if (isNChamberHit[0] && isNChamberHit[1]) {
        registry.fill(HIST("selected/hNhitsPerChamber_0"), 1., xCh[0], yCh[0]);
        registry.fill(HIST("selected/hNhitsPerChamber_1"), 1., xCh[1], yCh[1]);
      }
      if (isNChamberHit[0] && (!isNChamberHit[1])) {
        registry.fill(HIST("selected/hNhitsPerChamber_0"), 2., xCh[0], yCh[0]);
        registry.fill(HIST("selected/hNhitsPerChamber_1"), 3., xCh[1], yCh[1]);
      }
      if ((!isNChamberHit[0]) && isNChamberHit[1]) {
        registry.fill(HIST("selected/hNhitsPerChamber_0"), 3., xCh[0], yCh[0]);
        registry.fill(HIST("selected/hNhitsPerChamber_1"), 2., xCh[1], yCh[1]);
      }
      // Fill histos station 2, Chambers 3 & 4
      if (isNChamberHit[2] && isNChamberHit[3]) {
        registry.fill(HIST("selected/hNhitsPerChamber_2"), 1., xCh[2], yCh[2]);
        registry.fill(HIST("selected/hNhitsPerChamber_3"), 1., xCh[3], yCh[3]);
      }
      if (isNChamberHit[2] && (!isNChamberHit[3])) {
        registry.fill(HIST("selected/hNhitsPerChamber_2"), 2., xCh[2], yCh[2]);
        registry.fill(HIST("selected/hNhitsPerChamber_3"), 3., xCh[3], yCh[3]);
      }
      if ((!isNChamberHit[2]) && isNChamberHit[3]) {
        registry.fill(HIST("selected/hNhitsPerChamber_2"), 3., xCh[2], yCh[2]);
        registry.fill(HIST("selected/hNhitsPerChamber_3"), 2., xCh[3], yCh[3]);
      }
      // Fill histos station 3, Chambers 5 & 6
      if (isNChamberHit[4] && isNChamberHit[5]) {
        registry.fill(HIST("selected/hNhitsPerChamber_4"), 1., xCh[4], yCh[4]);
        registry.fill(HIST("selected/hNhitsPerChamber_5"), 1., xCh[5], yCh[5]);
      }
      if (isNChamberHit[4] && (!isNChamberHit[5])) {
        registry.fill(HIST("selected/hNhitsPerChamber_4"), 2., xCh[4], yCh[4]);
        registry.fill(HIST("selected/hNhitsPerChamber_5"), 3., xCh[5], yCh[5]);
      }
      if ((!isNChamberHit[4]) && isNChamberHit[5]) {
        registry.fill(HIST("selected/hNhitsPerChamber_4"), 3., xCh[4], yCh[4]);
        registry.fill(HIST("selected/hNhitsPerChamber_5"), 2., xCh[5], yCh[5]);
      }
      // Fill histos station 4, Chambers 7 & 8
      if (isNChamberHit[6] && isNChamberHit[7]) {
        registry.fill(HIST("selected/hNhitsPerChamber_6"), 1., xCh[6], yCh[6]);
        registry.fill(HIST("selected/hNhitsPerChamber_7"), 1., xCh[7], yCh[7]);
      }
      if (isNChamberHit[6] && (!isNChamberHit[7])) {
        registry.fill(HIST("selected/hNhitsPerChamber_6"), 2., xCh[6], yCh[6]);
        registry.fill(HIST("selected/hNhitsPerChamber_7"), 3., xCh[7], yCh[7]);
      }
      if ((!isNChamberHit[6]) && isNChamberHit[7]) {
        registry.fill(HIST("selected/hNhitsPerChamber_6"), 3., xCh[6], yCh[6]);
        registry.fill(HIST("selected/hNhitsPerChamber_7"), 2., xCh[7], yCh[7]);
      }
      // Fill histos station 5, Chambers 9 & 10
      if (isNChamberHit[8] && isNChamberHit[9]) {
        registry.fill(HIST("selected/hNhitsPerChamber_8"), 1., xCh[8], yCh[8]);
        registry.fill(HIST("selected/hNhitsPerChamber_9"), 1., xCh[9], yCh[9]);
      }
      if (isNChamberHit[8] && (!isNChamberHit[9])) {
        registry.fill(HIST("selected/hNhitsPerChamber_8"), 2., xCh[8], yCh[8]);
        registry.fill(HIST("selected/hNhitsPerChamber_9"), 3., xCh[9], yCh[9]);
      }
      if ((!isNChamberHit[8]) && isNChamberHit[9]) {
        registry.fill(HIST("selected/hNhitsPerChamber_8"), 3., xCh[8], yCh[8]);
        registry.fill(HIST("selected/hNhitsPerChamber_9"), 2., xCh[9], yCh[9]);
      }

    } // end filling info for selected candidates
  }

  /// process to fill histograms
  void FillHistosMC(double mEta, double mPhi, double mPt, uint16_t /*mchBitmap*/, bool isSel, double mGenEta, double mGenPt, double mGenPhi)
  {

    registry.fill(HIST("hPtRecPtGen"), mPt, mGenPt);
    registry.fill(HIST("hEtaRecEtaGen"), mEta, mGenEta);
    registry.fill(HIST("hPhiRecPhiGen"), mPhi, mGenPhi);

    /// fill histograms only for selected candidates
    ///
    if (isSel) {
      registry.fill(HIST("selected/hPtRecPtGen"), mPt, mGenPt);
      registry.fill(HIST("selected/hEtaRecEtaGen"), mEta, mGenEta);
      registry.fill(HIST("selected/hPhiRecPhiGen"), mPhi, mGenPhi);
    } // end filling info for selected candidates
  }

  /// Kinematic selection
  bool IsInKinematics(double eta, double pt)
  {
    bool isSelected = true;

    if (pt < ptMuonMin) {
      isSelected = false;
    }
    if ((eta < etaMuonMin) || (eta > etaMuonMax)) {
      isSelected = false;
    }
    return isSelected;
  }

  /// Event selection
  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent event)
  {
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(event); // extract event information and place it in the fValues array

    /// Analyse only selected events
    registry.fill(HIST("hEventCounter"), 0);
    if (!fEventCut->IsSelected(VarManager::fgValues)) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 1);
  }

  /// Muon selection and info filling
  template <uint32_t TMuonFillMap, typename TMuons>
  void runMuonSelection(TMuons const& tracksMuon)
  {
    /// loop on all muons
    LOGF(debug, " muon fwd tracks %i", tracksMuon.size());
    const int ncuts = fMuonCuts.size();
    std::vector<int> nselmuons(ncuts);
    for (int i = 0; i < ncuts; i++)
      nselmuons.push_back(0);

    rowCandidateBase.reserve(tracksMuon.size());
    for (auto& muon : tracksMuon) {

      VarManager::FillTrack<TMuonFillMap>(muon);

      LOGF(debug, "  %i / %f / %f / %f", muon.trackType(), muon.eta(), muon.pt(), muon.p());
      int mType = muon.trackType();
      double mPt = muon.pt();
      double mEta = TMath::Abs(muon.eta());
      double mPhi = muon.phi();
      uint16_t mchBitmap = muon.mchBitMap();

      /// select muons passing criteria
      std::vector<bool> isMuonSelected(ncuts);
      bool isMuonSelectedAny = false;
      int j = 0;

      ///   check the cuts and filters
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, j++) {
        isMuonSelected.push_back(false);
        LOGF(debug, " checking muon selected for cut %i", j);
        if ((*cut).IsSelected(VarManager::fgValues)) {
          LOGF(debug, " muon IS selected for cut %i", j);
          isMuonSelected[j] = true;
          isMuonSelectedAny = true;
          nselmuons[j] += 1;
        }
      }

      bool isKineAcc = IsInKinematics(mEta, mPt);
      if (!isKineAcc)
        isMuonSelectedAny = false;
      if (muonSelType >= 0) {
        if (mType != muonSelType)
          isMuonSelectedAny = false;
      }
      FillHistos(mEta, mPhi, mPt, mchBitmap, isMuonSelectedAny);

      if (isMuonSelectedAny)
        rowCandidateBase(mEta, mPt, mPhi, mchBitmap, mType);

    } // end loop on muons
    LOGF(debug, "end muon loop");

    // fill counter histo with number of muons with more than x muons selected
    if (nselmuons[0] > 0)
      registry.fill(HIST("hEventCounter"), 2);
    if (nselmuons[0] > 1)
      registry.fill(HIST("hEventCounter"), 3);
    if (nselmuons[0] > 2)
      registry.fill(HIST("hEventCounter"), 4);
    if (nselmuons[0] > 3)
      registry.fill(HIST("hEventCounter"), 5);
  }

  /// Muon selection and info filling
  template <uint32_t TMuonFillMap, typename TMuons>
  void runSimulatedMuonSelection(TMuons const& tracksMuon)
  {
    /// loop on all muons
    LOGF(debug, " muon fwd tracks %i", tracksMuon.size());
    const int ncuts = fMuonCuts.size();
    std::vector<int> nselmuons(ncuts);
    for (int i = 0; i < ncuts; i++)
      nselmuons.push_back(0);

    rowCandidateBase.reserve(tracksMuon.size());
    rowCandidateGen.reserve(tracksMuon.size());
    for (auto& muon : tracksMuon) {
      ///
      /// First compute MC matched quantities using either the DQ skimmed or the Framework data models
      double mGenPt = 0., mGenEta = 0., mGenPhi = 0.;
      LOGF(debug, " Looking for the correspondent MC particle");
      if constexpr ((TMuonFillMap & VarManager::ObjTypes::ReducedMuon) > 0) {
        auto particle = muon.reducedMCTrack();
        VarManager::FillTrack<gkParticleMCFillMap>(particle);

        auto pdgParticle = particle.pdgCode();
        if (!pdgParticle) {
          LOGF(warning, "MC particle PDG code not found, skip...");
          continue;
        }
        if (TMath::Abs(particle.pdgCode()) != muonPDG) {
          LOGF(warning, "MC particle does not correspond to a muon, skip...");
          continue;
        }
        mGenPt = particle.pt();
        mGenEta = TMath::Abs(particle.eta());
        mGenPhi = particle.phi();
      }
      if constexpr ((TMuonFillMap & VarManager::ObjTypes::Muon) > 0) {
        /// Check if the correspondent MC particle is really a muon
        ///   otherwise reject
        if (!muon.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        auto particle = muon.mcParticle();
        VarManager::FillTrack<gkParticleMCFillMap>(particle);
        LOGF(debug, " filled mc particle map");

        auto particlePdgCode = particle.pdgCode();

        if (TMath::Abs(particlePdgCode) != muonPDG) {
          LOGF(warning, "MC particle does not correspond to a muon, skip...");
          continue;
        }
        mGenPt = particle.pt();
        mGenEta = TMath::Abs(particle.eta());
        mGenPhi = particle.phi();
      }

      LOGF(debug, "now fill muon map");
      /// look the reconstructed quantities
      VarManager::FillTrack<TMuonFillMap>(muon);

      LOGF(debug, "  %i / %f / %f / %f", muon.trackType(), muon.eta(), muon.pt(), muon.p());
      int mType = muon.trackType();
      double mPt = muon.pt();
      double mEta = TMath::Abs(muon.eta());
      double mPhi = muon.phi();
      uint16_t mchBitmap = muon.mchBitMap();

      /// select muons passing criteria
      std::vector<bool> isMuonSelected(ncuts);
      bool isMuonSelectedAny = false;
      int j = 0;

      ///   check the cuts and filters
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, j++) {
        isMuonSelected.push_back(false);
        if ((*cut).IsSelected(VarManager::fgValues)) {
          isMuonSelected[j] = true;
          isMuonSelectedAny = true;
          nselmuons[j] += 1;
        }
      }

      bool isKineAcc = IsInKinematics(mEta, mPt);
      if (!isKineAcc)
        isMuonSelectedAny = false;
      if (muonSelType >= 0) {
        if (mType != muonSelType)
          isMuonSelectedAny = false;
      }
      /// fill histograms
      FillHistos(mEta, mPhi, mPt, mchBitmap, isMuonSelectedAny);
      FillHistosMC(mEta, mPhi, mPt, mchBitmap, isMuonSelectedAny, mGenEta, mGenPt, mGenPhi);

      if (isMuonSelectedAny) {
        rowCandidateBase(mEta, mPt, mPhi, mchBitmap, mType);
        rowCandidateGen(mGenEta, mGenPt, mGenPhi);
      }

    } // end loop on muons
    LOGF(debug, "end muon loop");

    // fill counter histo with number of muons with more than x muons selected
    if (nselmuons[0] > 0)
      registry.fill(HIST("hEventCounter"), 2);
    if (nselmuons[0] > 1)
      registry.fill(HIST("hEventCounter"), 3);
    if (nselmuons[0] > 2)
      registry.fill(HIST("hEventCounter"), 4);
    if (nselmuons[0] > 3)
      registry.fill(HIST("hEventCounter"), 5);
  }

  //! process function for full muon information
  void processReco(myEvents::iterator const& collision, aod::BCsWithTimestamps const&, myMuons const& muons)
  {
    /// Run event selection
    runEventSelection<gkEventFillMap>(collision);
    /// Run muon selection and fill output histograms
    runMuonSelection<gkMuonFillMap>(muons);
  }
  PROCESS_SWITCH(tableMakerMuonMchTrkEfficiency, processReco, "process reconstructed information", false);

  //! process function for reduced muon information
  void processRecoReduced(myReducedEvents::iterator const& event, myReducedMuons const& muons)
  {
    /// Run event selection
    runEventSelection<gkReducedEventFillMap>(event);
    /// Run muon selection and fill output histograms
    runMuonSelection<gkReducedMuonFillMap>(muons);
  }
  PROCESS_SWITCH(tableMakerMuonMchTrkEfficiency, processRecoReduced, "process reconstructed reduced information", true);

  //! process function for simulated muon information
  //! group according to reconstructed Collisions
  void processSim(myEventsMC::iterator const& collision, aod::BCsWithTimestamps const&, myMuonsMC const& muons,
                  aod::McParticles_001 const&, aod::McCollisions const&)
  {
    // TODO: investigate the collisions without corresponding mcCollision
    if (!collision.has_mcCollision()) {
      return;
    }
    /// Run event selection
    ///
    /// fill mc collision map
    auto mcCollision = collision.mcCollision();
    VarManager::FillEvent<gkEventMCFillMap>(mcCollision);
    ///  event selection
    runEventSelection<gkEventFillMap>(collision);

    /// Run muon selection and histo filling
    runSimulatedMuonSelection<gkMuonFillMap>(muons);
  }
  PROCESS_SWITCH(tableMakerMuonMchTrkEfficiency, processSim, "process simulation information", false);

  //! process function for reducedsimulated muon information
  //! group according to reconstructed Collisions
  void processSimReduced(myReducedEventsMC::iterator const& collision, myReducedMuonsMC const& muons,
                         aod::McParticles_001 const&, aod::McCollisions const&)
  {

    /// Run event selection
    ///
    /// fill mc collision map
    auto mcCollision = collision.reducedMCevent();
    VarManager::FillEvent<gkReducedEventMCFillMap>(mcCollision);
    /// event selection
    runEventSelection<gkReducedEventFillMap>(collision);

    /// Run muon selection and histo filling
    //        VarManager::ResetValues(0, VarManager::kNMCParticleVariables);
    runSimulatedMuonSelection<gkReducedMuonFillMap>(muons);
  }
  PROCESS_SWITCH(tableMakerMuonMchTrkEfficiency, processSimReduced, "process reconstructed reduced information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tableMakerMuonMchTrkEfficiency>(cfgc)};
}
