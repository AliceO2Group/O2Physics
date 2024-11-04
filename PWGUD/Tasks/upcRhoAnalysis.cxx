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
/// \brief  task for analysis of rho in UPCs using UD tables (from SG producer)
///         includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for decay phi anisotropy studies
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "random"
#include "TLorentzVector.h"

#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDSgCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

namespace o2::aod
{
namespace dipi
{
// general
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(NClass, nClass, int);
DECLARE_SOA_COLUMN(TotCharge, charge, int);
DECLARE_SOA_COLUMN(Pt, pT, double);
// system
DECLARE_SOA_COLUMN(M, m, double);
DECLARE_SOA_COLUMN(Rap, y, double);
DECLARE_SOA_COLUMN(PhiRandom, phiRandom, double);
DECLARE_SOA_COLUMN(PhiCharge, phiCharge, double);
// tracks
DECLARE_SOA_COLUMN(UdCollisionId, udCollisionId, int32_t);
DECLARE_SOA_COLUMN(Eta, eta, double);
DECLARE_SOA_COLUMN(Phi, phi, double);
DECLARE_SOA_COLUMN(Sign, sign, int);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, double);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, double);
DECLARE_SOA_COLUMN(NSigmaPi, nSigmaPi, double);
DECLARE_SOA_COLUMN(NSigmaEl, nSigmaEl, double);
} // namespace dipi
DECLARE_SOA_TABLE(SystemTree, "AOD", "SYSTEMTREE", dipi::RunNumber, dipi::NClass, dipi::TotCharge, dipi::M, dipi::Pt, dipi::Rap, dipi::PhiRandom, dipi::PhiCharge);
DECLARE_SOA_TABLE(TrackTree, "AOD", "TRACKTREE", dipi::RunNumber, dipi::NClass, dipi::UdCollisionId, dipi::Pt, dipi::Eta, dipi::Sign, dipi::DcaZ, dipi::DcaXY, dipi::NSigmaPi, dipi::NSigmaEl);
} // namespace o2::aod

struct upcRhoAnalysis {
  Produces<o2::aod::SystemTree> systemTree;
  Produces<o2::aod::TrackTree> trackTree;

  double PcEtaCut = 0.9; // physics coordination recommendation

  Configurable<bool> specifyGapSide{"specifyGapSide", true, "specify gap side for SG/DG produced data"};
  Configurable<int> gapSide{"gapSide", 2, "gap side for SG produced data"};
  Configurable<bool> requireTof{"requireTof", false, "require TOF signal"};

  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "max Z position cut on collisions"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZN common energy cut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZN time cut"};

  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "TPC nSigma pion cut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "max DCA cut on tracks"};
  Configurable<int> tracksMinItsNClsCut{"tracksMinItsNClsCut", 6, "min ITS clusters cut"};
  Configurable<double> tracksMaxItsChi2NClCut{"tracksMaxItsChi2NClCut", 3.0, "max ITS chi2/Ncls cut"};
  Configurable<int> tracksMinTpcNClsCut{"tracksMinTpcNClsCut", 120, "min TPC clusters cut"};
  Configurable<int> tracksMinTpcNClsCrossedRowsCut{"tracksMinTpcNClsCrossedRowsCut", 140, "min TPC crossed rows cut"};
  Configurable<double> tracksMinTpcChi2NClCut{"tracksMinTpcChi2NClCut", 1.0, "min TPC chi2/Ncls cut"};
  Configurable<double> tracksMaxTpcChi2NClCut{"tracksMaxTpcChi2NClCut", 1.8, "max TPC chi2/Ncls cut"};
  Configurable<double> tracksMinTpcNClsCrossedOverFindableCut{"tracksMinTpcNClsCrossedOverFindableCut", 1.05, "min TPC crossed rows / findable clusters cut"};
  Configurable<double> tracksMinPtCut{"tracksMinPtCut", 0.2, "min pT cut on tracks"};

  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.4, "min M cut for reco system"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 1.2, "max M cut for reco system"};
  Configurable<double> systemPtCut{"systemPtMaxCut", 0.1, "max pT cut for reco system"};
  Configurable<double> systemYCut{"systemYCut", 0.9, "rapiditiy cut for reco system"};

  ConfigurableAxis mAxis{"mAxis", {1000, 0.0, 10.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis mCutAxis{"mCutAxis", {160, 0.4, 1.2}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 10.0}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis ptCutAxis{"ptCutAxis", {100, 0.0, 0.1}, "p_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {100, 0.0, 0.01}, "p_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, o2::constants::math::TwoPI}, "#phi"};
  ConfigurableAxis phiAsymmAxis{"phiAsymmAxis", {182, -o2::constants::math::PI, o2::constants::math::PI}, "#phi"};
  ConfigurableAxis momentumFromPhiAxis{"momentumFromPhiAxis", {400, -0.1, 0.1}, "p (GeV/#it{c})"};
  ConfigurableAxis ptQuantileAxis{"ptQuantileAxis", {0, 0.0181689, 0.0263408, 0.0330488, 0.0390369, 0.045058, 0.0512604, 0.0582598, 0.066986, 0.0788085, 0.1}, "p_{T} (GeV/#it{c})"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // QA //
    // collisions
    registry.add("QC/collisions/hPosXY", ";x (cm);y (cm);counts", kTH2D, {{2000, -0.1, 0.1}, {2000, -0.1, 0.1}});
    registry.add("QC/collisions/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/collisions/hNumContrib", ";number of contributors;counts", kTH1D, {{36, -0.5, 35.5}});
    registry.add("QC/collisions/hZdcCommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0}, {250, -5.0, 20.0}});
    registry.add("QC/collisions/hZdcTime", ";ZNA time (ns);ZNC time (ns);counts", kTH2D, {{200, -10.0, 10.0}, {200, -10.0, 10.0}});
    // all tracks
    registry.add("QC/tracks/raw/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {{400, -10.0, 30.0}});
    registry.add("QC/tracks/raw/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/tracks/raw/hTpcNSigmaEl", ";TPC n#sigma_{e};counts", kTH1D, {{400, -10.0, 30.0}});
    registry.add("QC/tracks/raw/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    registry.add("QC/tracks/raw/hItsNCls", ";ITS N_{cls};counts", kTH1D, {{11, -0.5, 10.5}});
    registry.add("QC/tracks/raw/hItsChi2NCl", ";ITS #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    registry.add("QC/tracks/raw/hTpcChi2NCl", ";TPC #chi^{2}/N_{cls};counts", kTH1D, {{1000, 0.0, 100.0}});
    registry.add("QC/tracks/raw/hTpcNCls", ";TPC N_{cls} findable;counts", kTH1D, {{200, 0.0, 200.0}});
    registry.add("QC/tracks/raw/hTpcNClsCrossedRows", ";TPC crossed rows;counts", kTH1D, {{200, 0.0, 200.0}});
    // track quality selections vs system mass
    registry.add("QC/tracks/2D/mass/leading/hItsNClsVsM", ";m (GeV/#it{c}^{2});ITS N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/mass/leading/hItsChi2NClVsM", ";m (GeV/#it{c}^{2});ITS #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/mass/leading/hTpcChi2NClVsM", ";m (GeV/#it{c}^{2});TPC #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/mass/leading/hTpcNClsVsM", ";m (GeV/#it{c}^{2});TPC N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/mass/leading/hTpcNClsCrossedRowsVsM", ";m (GeV/#it{c}^{2});TPC crossed rows of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/mass/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsM", ";m (GeV/#it{c}^{2});TPC crossed rows / TPC N_{cls} findable of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 10.0}});
    registry.add("QC/tracks/2D/mass/subleading/hItsNClsVsM", ";m (GeV/#it{c}^{2});ITS N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/mass/subleading/hItsChi2NClVsM", ";m (GeV/#it{c}^{2});ITS #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/mass/subleading/hTpcChi2NClVsM", ";m (GeV/#it{c}^{2});TPC #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/mass/subleading/hTpcNClsVsM", ";m (GeV/#it{c}^{2});TPC N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/mass/subleading/hTpcNClsCrossedRowsVsM", ";m (GeV/#it{c}^{2});TPC crossed rows of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/mass/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsM", ";m (GeV/#it{c}^{2});TPC crossed rows / TPC N_{cls} findable of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 10.0}});
    // track quality selections vs system rapidity
    registry.add("QC/tracks/2D/rapidity/leading/hItsNClsVsY", ";y;ITS N_{cls} of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/rapidity/leading/hItsChi2NClVsY", ";y;ITS #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/rapidity/leading/hTpcChi2NClVsY", ";y;TPC #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/rapidity/leading/hTpcNClsVsY", ";y;TPC N_{cls} of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/rapidity/leading/hTpcNClsCrossedRowsVsY", ";y;TPC crossed rows of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/rapidity/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsY", ";y;TPC crossed rows / TPC N_{cls} findable of leading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 10.0}});
    registry.add("QC/tracks/2D/rapidity/subleading/hItsNClsVsY", ";y;ITS N_{cls} of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/rapidity/subleading/hItsChi2NClVsY", ";y;ITS #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/rapidity/subleading/hTpcChi2NClVsY", ";y;TPC #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/rapidity/subleading/hTpcNClsVsY", ";y;TPC N_{cls} of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/rapidity/subleading/hTpcNClsCrossedRowsVsY", ";y;TPC crossed rows of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/rapidity/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsY", ";y;TPC crossed rows / TPC N_{cls} findable of subleading-#it{p} track;counts", kTH2D, {{180, -0.9, 0.9}, {1000, 0.0, 10.0}});
    // track quality selections vs system pT
    registry.add("QC/tracks/2D/pT/leading/hItsNClsVsPt", ";p_{T} (GeV/#it{c});ITS N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/pT/leading/hItsChi2NClVsPt", ";p_{T} (GeV/#it{c});ITS #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/pT/leading/hTpcChi2NClVsPt", ";p_{T} (GeV/#it{c});TPC #chi^{2}/N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/pT/leading/hTpcNClsVsPt", ";p_{T} (GeV/#it{c});TPC N_{cls} of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/pT/leading/hTpcNClsCrossedRowsVsPt", ";p_{T} (GeV/#it{c});TPC crossed rows of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/pT/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsPt", ";p_{T} (GeV/#it{c});TPC crossed rows / TPC N_{cls} findable of leading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 10.0}});
    registry.add("QC/tracks/2D/pT/subleading/hItsNClsVsPt", ";p_{T} (GeV/#it{c});ITS N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {11, -0.5, 10.5}});
    registry.add("QC/tracks/2D/pT/subleading/hItsChi2NClVsPt", ";p_{T} (GeV/#it{c});ITS #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/pT/subleading/hTpcChi2NClVsPt", ";p_{T} (GeV/#it{c});TPC #chi^{2}/N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 100.0}});
    registry.add("QC/tracks/2D/pT/subleading/hTpcNClsVsPt", ";p_{T} (GeV/#it{c});TPC N_{cls} of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/pT/subleading/hTpcNClsCrossedRowsVsPt", ";p_{T} (GeV/#it{c});TPC crossed rows of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {200, 0.0, 200.0}});
    registry.add("QC/tracks/2D/pT/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsPt", ";p_{T} (GeV/#it{c});TPC crossed rows / TPC N_{cls} findable of subleading-#it{p} track;counts", kTH2D, {{1000, 0.0, 10.0}, {1000, 0.0, 10.0}});
    // tracks passing selections
    registry.add("QC/tracks/cut/hTpcNSigmaPi2D", ";TPC n#sigma(#pi_{leading});TPC n#sigma(#pi_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    registry.add("QC/tracks/cut/hTpcNSigmaEl2D", ";TPC n#sigma(e_{leading});TPC n#sigma(e_{subleading});counts", kTH2D, {{400, -10.0, 30.0}, {400, -10.0, 30.0}});
    registry.add("QC/tracks/cut/hTpcSignalVsP", ";p (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    registry.add("QC/tracks/cut/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    registry.add("QC/tracks/cut/hRemainingTracks", ";remaining tracks;counts", kTH1D, {{21, -0.5, 20.5}});
    registry.add("QC/tracks/cut/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {{1000, -5.0, 5.0}, {1000, -5.0, 5.0}});
    // selection counter
    std::vector<std::string> selectionCounterLabels = {"all tracks", "PV contributor", "ITS hit", "ITS N_{clusters}", "ITS #chi^{2}/N_{clusters}", "TPC hit", "TPC N_{clusters}", "TPC #chi^{2}/N_{clusters}", "TPC crossed rows",
                                                       "TPC crossed rows/N_{clusters}"
                                                       "TOF requirement",
                                                       "p_{T}", "DCA", "#eta", "exactly 2 tracks", "PID"};
    auto hSelectionCounter = registry.add<TH1>("QC/tracks/hSelectionCounter", ";;counts", kTH1D, {{static_cast<int>(selectionCounterLabels.size()), -0.5, static_cast<double>(selectionCounterLabels.size()) - 0.5}});
    for (int i = 0; i < static_cast<int>(selectionCounterLabels.size()); ++i)
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, selectionCounterLabels[i].c_str());
    // TOF hit check
    registry.add("QC/tracks/hTofHitCheck", ";leading track TOF hit;subleading track TOF hit;counts", kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    // RECO HISTOS //
    // PIONS
    // no selection
    registry.add("pions/no-selection/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("pions/no-selection/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("pions/no-selection/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("pions/no-selection/like-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("pions/no-selection/like-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("pions/no-selection/like-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    // selected
    registry.add("pions/selected/unlike-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("pions/selected/unlike-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("pions/selected/unlike-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("pions/selected/like-sign/hPt", ";p_{T}(#pi_{leading}) (GeV/#it{c});p_{T}(#pi_{subleading}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("pions/selected/like-sign/hEta", ";#eta(#pi_{leading});#eta(#pi_{subleading});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("pions/selected/like-sign/hPhi", ";#phi(#pi_{leading});#phi(#pi_{subleading});counts", kTH2D, {phiAxis, phiAxis});

    // RAW RHOS
    registry.add("system/raw/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("system/raw/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("system/raw/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("system/raw/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/raw/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("system/raw/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("system/raw/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("system/raw/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/raw/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("system/raw/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("system/raw/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("system/raw/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});

    // SELECTED RHOS
    // no selection
    registry.add("system/cut/no-selection/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/no-selection/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/no-selection/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/unlike-sign/hTofClassVsM", ";m (GeV/#it{c}^{2});TOF class;counts", kTH2D, {mCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/unlike-sign/hTofClassVsPt", ";p_{T} (GeV/#it{c});TOF class;counts", kTH2D, {ptCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/unlike-sign/hTofClassVsY", ";y;TOF class;counts", kTH2D, {yAxis, {4, -0.5, 3.5}});

    registry.add("system/cut/no-selection/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/no-selection/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/like-sign/positive/hTofClassVsM", ";m (GeV/#it{c}^{2});TOF class;counts", kTH2D, {mCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/like-sign/positive/hTofClassVsPt", ";p_{T} (GeV/#it{c});TOF class;counts", kTH2D, {ptCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/like-sign/positive/hTofClassVsY", ";y;TOF class;counts", kTH2D, {yAxis, {4, -0.5, 3.5}});

    registry.add("system/cut/no-selection/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/no-selection/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/no-selection/like-sign/negative/hTofClassVsM", ";m (GeV/#it{c}^{2});TOF class;counts", kTH2D, {mCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/like-sign/negative/hTofClassVsPt", ";p_{T} (GeV/#it{c});TOF class;counts", kTH2D, {ptCutAxis, {4, -0.5, 3.5}});
    registry.add("system/cut/no-selection/like-sign/negative/hTofClassVsY", ";y;TOF class;counts", kTH2D, {yAxis, {4, -0.5, 3.5}});

    // 0n0n
    registry.add("system/cut/0n0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0n0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0n0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0n0n/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/0n0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0n0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0n0n/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/0n0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0n0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0n0n/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // Xn0n
    registry.add("system/cut/Xn0n/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/Xn0n/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/Xn0n/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/Xn0n/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/Xn0n/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/Xn0n/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/Xn0n/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // 0nXn
    registry.add("system/cut/0nXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0nXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0nXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0nXn/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/0nXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0nXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0nXn/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/0nXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/0nXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/0nXn/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    // XnXn
    registry.add("system/cut/XnXn/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/XnXn/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/XnXn/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/XnXn/unlike-sign/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/XnXn/like-sign/positive/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/XnXn/like-sign/positive/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/XnXn/like-sign/positive/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});

    registry.add("system/cut/XnXn/like-sign/negative/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mCutAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptCutAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPt2", ";p_{T}^{2} (GeV^{2}/#it{c}^{2});counts", kTH1D, {pt2Axis});
    registry.add("system/cut/XnXn/like-sign/negative/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mCutAxis, ptCutAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPhiRandom", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPhiCharge", ";#phi;counts", kTH1D, {phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPhiRandomVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPhiChargeVsM", ";m (GeV/#it{c}^{2});#phi;counts", kTH2D, {mCutAxis, phiAsymmAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPyVsPxRandom", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
    registry.add("system/cut/XnXn/like-sign/negative/hPyVsPxCharge", ";p_{x} (GeV/#it{c});p_{y} (GeV/#it{c});counts", kTH2D, {momentumFromPhiAxis, momentumFromPhiAxis});
  }

  template <typename C>
  bool collisionPassesCuts(const C& collision) // collision cuts
  {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut)
      return false;
    if (specifyGapSide && collision.gapSide() != gapSide)
      return false;
    return true;
  }

  template <typename T>
  bool trackPassesCuts(const T& track) // track cuts (PID done separately)
  {
    if (!track.isPVContributor())
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 1);

    if (!track.hasITS())
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 2);

    if (track.itsNCls() < tracksMinItsNClsCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 3);

    if (track.itsChi2NCl() > tracksMaxItsChi2NClCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 4);

    if (!track.hasTPC())
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 5);

    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < tracksMinTpcNClsCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 6);

    if (track.tpcChi2NCl() > tracksMaxTpcChi2NClCut || track.tpcChi2NCl() < tracksMinTpcChi2NClCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 7);

    if (track.tpcNClsCrossedRows() < tracksMinTpcNClsCrossedRowsCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 8);

    if (static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()) < tracksMinTpcNClsCrossedOverFindableCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 9);

    if (requireTof && !track.hasTOF())
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 10);

    if (track.pt() < tracksMinPtCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 11);

    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > (0.0105 + 0.0350 / std::pow(track.pt(), 1.01)))
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 12);

    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut)
      return false;
    registry.fill(HIST("QC/tracks/hSelectionCounter"), 13);
    // if all selections passed
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) // n-dimensional PID cut
  {
    double radius = 0.0;
    for (const auto& track : cutTracks)
      radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  int tracksTotalCharge(const T& cutTracks) // total charge of selected tracks
  {
    int charge = 0;
    for (const auto& track : cutTracks)
      charge += track.sign();
    return charge;
  }

  bool systemPassCuts(const TLorentzVector& system) // system cuts
  {
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut)
      return false;
    if (system.Pt() > systemPtCut)
      return false;
    if (std::abs(system.Rapidity()) > systemYCut)
      return false;
    return true;
  }

  TLorentzVector reconstructSystem(const std::vector<TLorentzVector>& cutTracks4Vecs) // reconstruct system from 4-vectors
  {
    TLorentzVector system;
    for (const auto& track4Vec : cutTracks4Vecs)
      system += track4Vec;
    return system;
  }

  double getPhiRandom(const std::vector<TLorentzVector>& cutTracks4Vecs) // decay phi anisotropy
  {                                                                      // two possible definitions of phi: randomize the tracks
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    TLorentzVector pOne = cutTracks4Vecs[indices[0]];
    TLorentzVector pTwo = cutTracks4Vecs[indices[1]];
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  template <typename T>
  double getPhiCharge(const T& cutTracks, const std::vector<TLorentzVector>& cutTracks4Vecs)
  { // two possible definitions of phi: charge-based assignment
    TLorentzVector pOne, pTwo;
    if (cutTracks[0].sign() > 0) {
      pOne = cutTracks4Vecs[0];
      pTwo = cutTracks4Vecs[1];
    } else {
      pOne = cutTracks4Vecs[1];
      pTwo = cutTracks4Vecs[0];
    }
    TLorentzVector pPlus = pOne + pTwo;
    TLorentzVector pMinus = pOne - pTwo;
    return pPlus.DeltaPhi(pMinus);
  }

  void processReco(FullUDSgCollision const& collision, FullUDTracks const& tracks)
  {
    // QC histograms
    registry.fill(HIST("QC/collisions/hPosXY"), collision.posX(), collision.posY());
    registry.fill(HIST("QC/collisions/hPosZ"), collision.posZ());
    registry.fill(HIST("QC/collisions/hZdcCommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    registry.fill(HIST("QC/collisions/hZdcTime"), collision.timeZNA(), collision.timeZNC());
    registry.fill(HIST("QC/collisions/hNumContrib"), collision.numContrib());

    if (!collisionPassesCuts(collision))
      return;

    // event tagging
    bool XnXn = false, OnOn = false, XnOn = false, OnXn = false; // note: On == 0n...
    int nClass = -1;
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() < ZNcommonEnergyCut) {
      OnOn = true;
      nClass = 0;
    }
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut &&
        collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      XnXn = true;
      nClass = 3;
    }
    if (collision.energyCommonZNA() > ZNcommonEnergyCut && std::abs(collision.timeZNA()) < ZNtimeCut && collision.energyCommonZNC() < ZNcommonEnergyCut) {
      XnOn = true;
      nClass = 1;
    }
    if (collision.energyCommonZNA() < ZNcommonEnergyCut && collision.energyCommonZNC() > ZNcommonEnergyCut && std::abs(collision.timeZNC()) < ZNtimeCut) {
      OnXn = true;
      nClass = 2;
    }
    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<TLorentzVector> cutTracks4Vecs;

    int trackCounter = 0;
    for (const auto& track : tracks) {
      // double p = momentum(track.px(), track.py(), track.pz());
      registry.fill(HIST("QC/tracks/raw/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/tracks/raw/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/tracks/raw/hTpcNSigmaEl"), track.tpcNSigmaEl());
      registry.fill(HIST("QC/tracks/raw/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      registry.fill(HIST("QC/tracks/raw/hItsNCls"), track.itsNCls());
      registry.fill(HIST("QC/tracks/raw/hItsChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("QC/tracks/raw/hTpcChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("QC/tracks/raw/hTpcNCls"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      registry.fill(HIST("QC/tracks/raw/hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("QC/tracks/hSelectionCounter"), 0);

      if (!trackPassesCuts(track))
        continue;
      trackCounter++;
      cutTracks.push_back(track);
      TLorentzVector track4Vec;
      track4Vec.SetXYZM(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged); // apriori assume pion mass
      cutTracks4Vecs.push_back(track4Vec);
      registry.fill(HIST("QC/tracks/cut/hTpcSignalVsP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
      registry.fill(HIST("QC/tracks/cut/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("QC/tracks/cut/hDcaXYZ"), track.dcaZ(), track.dcaXY());
    }
    registry.fill(HIST("QC/tracks/cut/hRemainingTracks"), trackCounter);

    if (cutTracks.size() != cutTracks4Vecs.size()) {
      LOG(error);
      return;
    }

    if (cutTracks.size() != 2)
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      registry.fill(HIST("QC/tracks/hSelectionCounter"), 14);

    registry.fill(HIST("QC/tracks/cut/hTpcNSigmaPi2D"), cutTracks[0].tpcNSigmaPi(), cutTracks[1].tpcNSigmaPi());
    for (int i = 0; i <= 1; i++)
      trackTree(collision.runNumber(), nClass, cutTracks[i].udCollisionId(), cutTracks[i].pt(), eta(cutTracks[i].px(), cutTracks[i].py(), cutTracks[i].pz()), cutTracks[i].sign(), cutTracks[i].dcaZ(), cutTracks[i].dcaXY(), cutTracks[i].tpcNSigmaPi(), cutTracks[i].tpcNSigmaEl());

    if (!tracksPassPiPID(cutTracks))
      return;
    for (int i = 0; i < static_cast<int>(cutTracks.size()); i++)
      registry.fill(HIST("QC/tracks/hSelectionCounter"), 15);
    registry.fill(HIST("QC/tracks/cut/hTpcNSigmaEl2D"), cutTracks[0].tpcNSigmaEl(), cutTracks[1].tpcNSigmaEl());

    // reonstruct system and calculate total charge, save commonly used values into variables
    TLorentzVector system = reconstructSystem(cutTracks4Vecs);
    int totalCharge = tracksTotalCharge(cutTracks);
    double mass = system.M();
    double pT = system.Pt();
    double pTsquare = pT * pT;
    double rapidity = system.Rapidity();
    double phiRandom = getPhiRandom(cutTracks4Vecs);
    double phiCharge = getPhiCharge(cutTracks, cutTracks4Vecs);
    // differentiate leading- and subleading-momentum tracks
    auto leadingMomentumTrack = momentum(cutTracks[0].px(), cutTracks[0].py(), cutTracks[0].pz()) > momentum(cutTracks[1].px(), cutTracks[1].py(), cutTracks[1].pz()) ? cutTracks[0] : cutTracks[1];
    auto subleadingMomentumTrack = (leadingMomentumTrack == cutTracks[0]) ? cutTracks[1] : cutTracks[0];
    double leadingPt = leadingMomentumTrack.pt();
    double subleadingPt = subleadingMomentumTrack.pt();
    double leadingEta = eta(leadingMomentumTrack.px(), leadingMomentumTrack.py(), leadingMomentumTrack.pz());
    double subleadingEta = eta(subleadingMomentumTrack.px(), subleadingMomentumTrack.py(), subleadingMomentumTrack.pz());
    double leadingPhi = phi(leadingMomentumTrack.px(), leadingMomentumTrack.py());
    double subleadingPhi = phi(subleadingMomentumTrack.px(), subleadingMomentumTrack.py());
    // fill TOF hit checker
    registry.fill(HIST("QC/tracks/hTofHitCheck"), leadingMomentumTrack.hasTOF(), subleadingMomentumTrack.hasTOF());
    int tofClass = -1;
    if (!leadingMomentumTrack.hasTOF() && !subleadingMomentumTrack.hasTOF())
      tofClass = 0;
    else if (leadingMomentumTrack.hasTOF() && !subleadingMomentumTrack.hasTOF())
      tofClass = 1;
    else if (!leadingMomentumTrack.hasTOF() && subleadingMomentumTrack.hasTOF())
      tofClass = 2;
    else if (leadingMomentumTrack.hasTOF() && subleadingMomentumTrack.hasTOF())
      tofClass = 3;
    // fill 2D track QC histograms
    // mass
    registry.fill(HIST("QC/tracks/2D/mass/leading/hItsNClsVsM"), mass, leadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/mass/leading/hItsChi2NClVsM"), mass, leadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/mass/leading/hTpcChi2NClVsM"), mass, leadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/mass/leading/hTpcNClsVsM"), mass, leadingMomentumTrack.tpcNClsFindable() - leadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/mass/leading/hTpcNClsCrossedRowsVsM"), mass, leadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/mass/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsM"), mass, (static_cast<float>(leadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(leadingMomentumTrack.tpcNClsFindable())));
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hItsNClsVsM"), mass, subleadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hItsChi2NClVsM"), mass, subleadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hTpcChi2NClVsM"), mass, subleadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hTpcNClsVsM"), mass, subleadingMomentumTrack.tpcNClsFindable() - subleadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hTpcNClsCrossedRowsVsM"), mass, subleadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/mass/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsM"), mass, (static_cast<float>(subleadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(subleadingMomentumTrack.tpcNClsFindable())));
    // rapidity
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hItsNClsVsY"), rapidity, leadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hItsChi2NClVsY"), rapidity, leadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hTpcChi2NClVsY"), rapidity, leadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hTpcNClsVsY"), rapidity, leadingMomentumTrack.tpcNClsFindable() - leadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hTpcNClsCrossedRowsVsY"), rapidity, leadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/rapidity/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsY"), rapidity, (static_cast<float>(leadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(leadingMomentumTrack.tpcNClsFindable())));
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hItsNClsVsY"), rapidity, subleadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hItsChi2NClVsY"), rapidity, subleadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hTpcChi2NClVsY"), rapidity, subleadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hTpcNClsVsY"), rapidity, subleadingMomentumTrack.tpcNClsFindable() - subleadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hTpcNClsCrossedRowsVsY"), rapidity, subleadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/rapidity/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsY"), rapidity, (static_cast<float>(subleadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(subleadingMomentumTrack.tpcNClsFindable())));
    // pT
    registry.fill(HIST("QC/tracks/2D/pT/leading/hItsNClsVsPt"), pT, leadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/pT/leading/hItsChi2NClVsPt"), pT, leadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/pT/leading/hTpcChi2NClVsPt"), pT, leadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/pT/leading/hTpcNClsVsPt"), pT, leadingMomentumTrack.tpcNClsFindable() - leadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/pT/leading/hTpcNClsCrossedRowsVsPt"), pT, leadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/pT/leading/hTpcNClsCrossedRowsOverTpcNClsFindableVsPt"), pT, (static_cast<float>(leadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(leadingMomentumTrack.tpcNClsFindable())));
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hItsNClsVsPt"), pT, subleadingMomentumTrack.itsNCls());
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hItsChi2NClVsPt"), pT, subleadingMomentumTrack.itsChi2NCl());
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hTpcChi2NClVsPt"), pT, subleadingMomentumTrack.tpcChi2NCl());
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hTpcNClsVsPt"), pT, subleadingMomentumTrack.tpcNClsFindable() - subleadingMomentumTrack.tpcNClsFindableMinusFound());
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hTpcNClsCrossedRowsVsPt"), pT, subleadingMomentumTrack.tpcNClsCrossedRows());
    registry.fill(HIST("QC/tracks/2D/pT/subleading/hTpcNClsCrossedRowsOverTpcNClsFindableVsPt"), pT, (static_cast<float>(subleadingMomentumTrack.tpcNClsCrossedRows()) / static_cast<float>(subleadingMomentumTrack.tpcNClsFindable())));
    // fill tree
    systemTree(collision.runNumber(), nClass, totalCharge, mass, pT, rapidity, phiRandom, phiCharge);
    // fill raw histograms according to the total charge
    switch (totalCharge) {
      case 0:
        registry.fill(HIST("pions/no-selection/unlike-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/no-selection/unlike-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/no-selection/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/raw/unlike-sign/hM"), mass);
        registry.fill(HIST("system/raw/unlike-sign/hPt"), pT);
        registry.fill(HIST("system/raw/unlike-sign/hPtVsM"), mass, pT);
        registry.fill(HIST("system/raw/unlike-sign/hY"), rapidity);
        break;

      case 2:
        registry.fill(HIST("pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/raw/like-sign/positive/hM"), mass);
        registry.fill(HIST("system/raw/like-sign/positive/hPt"), pT);
        registry.fill(HIST("system/raw/like-sign/positive/hPtVsM"), mass, pT);
        registry.fill(HIST("system/raw/like-sign/positive/hY"), rapidity);
        break;

      case -2:
        registry.fill(HIST("pions/no-selection/like-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/no-selection/like-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/no-selection/like-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/raw/like-sign/negative/hM"), mass);
        registry.fill(HIST("system/raw/like-sign/negative/hPt"), pT);
        registry.fill(HIST("system/raw/like-sign/negative/hPtVsM"), mass, pT);
        registry.fill(HIST("system/raw/like-sign/negative/hY"), rapidity);
        break;

      default:
        break;
    }

    // apply cuts to system
    if (!systemPassCuts(system))
      return;

    // fill histograms for system passing cuts
    switch (totalCharge) {
      case 0:
        registry.fill(HIST("pions/selected/unlike-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/selected/unlike-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/selected/unlike-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hM"), mass);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPt"), pT);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPt2"), pTsquare);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPtVsM"), mass, pT);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hY"), rapidity);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandom"), phiRandom);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPhiCharge"), phiCharge);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hTofClassVsM"), mass, tofClass);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hTofClassVsPt"), pT, tofClass);
        registry.fill(HIST("system/cut/no-selection/unlike-sign/hTofClassVsY"), rapidity, tofClass);
        if (OnOn) {
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hM"), mass);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPt"), pT);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hY"), rapidity);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0n0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hM"), mass);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPt"), pT);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hY"), rapidity);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/Xn0n/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hM"), mass);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPt"), pT);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hY"), rapidity);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0nXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hM"), mass);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPt"), pT);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hY"), rapidity);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/XnXn/unlike-sign/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case 2:
        registry.fill(HIST("pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hM"), mass);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPt"), pT);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPt2"), pTsquare);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPtVsM"), mass, pT);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hY"), rapidity);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandom"), phiRandom);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiCharge"), phiCharge);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hTofClassVsM"), mass, tofClass);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hTofClassVsPt"), pT, tofClass);
        registry.fill(HIST("system/cut/no-selection/like-sign/positive/hTofClassVsY"), rapidity, tofClass);
        if (OnOn) {
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hM"), mass);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPt"), pT);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hY"), rapidity);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0n0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hM"), mass);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt"), pT);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hY"), rapidity);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/Xn0n/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hM"), mass);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPt"), pT);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hY"), rapidity);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0nXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hM"), mass);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPt"), pT);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hY"), rapidity);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/XnXn/like-sign/positive/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      case -2:
        registry.fill(HIST("pions/selected/like-sign/hPt"), leadingPt, subleadingPt);
        registry.fill(HIST("pions/selected/like-sign/hEta"), leadingEta, subleadingEta);
        registry.fill(HIST("pions/selected/like-sign/hPhi"), leadingPhi, subleadingPhi);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hM"), mass);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPt"), pT);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPt2"), pTsquare);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPtVsM"), mass, pT);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hY"), rapidity);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandom"), phiRandom);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiCharge"), phiCharge);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hTofClassVsM"), mass, tofClass);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hTofClassVsPt"), pT, tofClass);
        registry.fill(HIST("system/cut/no-selection/like-sign/negative/hTofClassVsY"), rapidity, tofClass);
        if (OnOn) {
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hM"), mass);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPt"), pT);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hY"), rapidity);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0n0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnOn) {
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hM"), mass);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt"), pT);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hY"), rapidity);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/Xn0n/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (OnXn) {
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hM"), mass);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPt"), pT);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hY"), rapidity);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/0nXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        } else if (XnXn) {
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hM"), mass);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPt"), pT);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPt2"), pTsquare);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPtVsM"), mass, pT);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hY"), rapidity);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandom"), phiRandom);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiCharge"), phiCharge);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiRandomVsM"), mass, phiRandom);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPhiChargeVsM"), mass, phiCharge);
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxRandom"), pT * std::cos(phiRandom), pT * std::sin(phiRandom));
          registry.fill(HIST("system/cut/XnXn/like-sign/negative/hPyVsPxCharge"), pT * std::cos(phiCharge), pT * std::sin(phiCharge));
        }
        break;

      default:
        break;
    }
  }
  PROCESS_SWITCH(upcRhoAnalysis, processReco, "analyse reco tracks", true);

  // void processMC(aod::UDMcCollisions::iterator const&, aod::UDMcParticles const& mcparticles)
  // {
  //   // loop over all particles in the event
  //   for (auto const& mcparticle : mcparticles) {
  //     // only consider charged pions
  //     if (std::abs(mcparticle.pdgCode()) != 211)
  //       continue;
  //   }
  // }
  // PROCESS_SWITCH(upcRhoAnalysis, processMC, "analyse MC tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)};
}
