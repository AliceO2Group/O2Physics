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

/// \file dataCreatorJpsiHadReduced.cxx
/// \brief Creation of J/Psi-LF hadron pairs for Beauty hadron analyses
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/WorkflowSpec.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

enum Event : uint8_t { // TODO: check if needed
  Processed = 0,
  NoCharmHadPiSelected,
  CharmHadPiSelected,
  NEvent
};

enum DecayChannel : uint8_t {
  B0ToJpsiK0Star = 0,
  BplusToJpsiK,
  BsToJpsiPhi
};

enum WrongCollisionType : uint8_t {
  None = 0,
  WrongAssociation,
  SplitCollision,
};

/// Creation of Jpsi-Had pairs for Beauty hadrons
struct HfDataCreatorJpsiHadReduced {
  // Produces AOD tables to store track information
  // collision related tables
  Produces<aod::HfRedCollisions> hfReducedCollision;
  Produces<aod::HfRedCollCents> hfReducedCollCentrality;
  Produces<aod::HfRedQvectors> hfReducedQvector;
  Produces<aod::HfRedCollExtras> hfReducedCollExtra;
  Produces<aod::HfOrigColCounts> hfCollisionCounter;
  // J/Psi related tables
  Produces<aod::HfRedJpsis> hfJpsi;
  Produces<aod::HfRedJpsiCov> hfRedJpsiCov;
  // Ka bachelor related tables
  Produces<aod::HfRedBach0Bases> hfTrackLfDau0;
  Produces<aod::HfRedBach0Cov> hfTrackCovLfDau0;
  Produces<aod::HfRedBach1Bases> hfTrackLfDau1;
  Produces<aod::HfRedBach1Cov> hfTrackCovLfDau1;
  // MC related tables
  Produces<aod::HfMcRecRedJPKs> rowHfJpsiKMcRecReduced;
  Produces<aod::HfMcRecRedJPPhis> rowHfJpsiPhiMcRecReduced;
  Produces<aod::HfMcGenRedBps> rowHfBpMcGenReduced;
  Produces<aod::HfMcGenRedBss> rowHfBsMcGenReduced;

  Produces<aod::HfCfgBpToJpsi> rowCandidateConfigBplus;
  Produces<aod::HfCfgBsToJpsis> rowCandidateConfigBs;

  Configurable<bool> skipRejectedCollisions{"skipRejectedCollisions", true, "skips collisions rejected by the event selection, instead of flagging only"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};

  Configurable<bool> runJpsiToee{"runJpsiToee", false, "Run analysis for J/Psi to ee (debug)"};
  struct : o2::framework::ConfigurableGroup {
    // TPC PID
    Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
    Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
    Configurable<double> nSigmaTpcElMax{"nSigmaTpcElMax", 4., "Electron nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcPiMin{"nSigmaTpcPiMin", 2.5, "Pion nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcPiMax{"nSigmaTpcPiMax", 99., "Pion nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcPrMin{"nSigmaTpcPrMin", -99., "Proton nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcPrMax{"nSigmaTpcPrMax", 99, "Proton nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcElCombinedMax{"nSigmaTpcElCombinedMax", 4., "Electron Nsigma cut on TPC combined with TOF"};
    Configurable<double> nSigmaTpcPiCombinedMin{"nSigmaTpcPiCombinedMin", 2.5, "Pion Nsigma cut on TPC combined with TOF"};
    Configurable<double> nSigmaTpcPiCombinedMax{"nSigmaTpcPiCombinedMax", 99., "Pion Nsigma cut on TPC combined with TOF"};
    Configurable<double> nSigmaTpcPrCombinedMin{"nSigmaTpcPrCombinedMin", -99., "Proton Nsigma cut on TPC combined with TOF"};
    Configurable<double> nSigmaTpcPrCombinedMax{"nSigmaTpcPrCombinedMax", 99, "Proton Nsigma cut on TPC combined with TOF"};
    // TOF PID
    Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
    Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
    Configurable<double> nSigmaTofElMax{"nSigmaTofElMax", 4., "Electron nsigma cut on TPC only"};
    Configurable<double> nSigmaTofPiMin{"nSigmaTofPiMin", -99, "Pion nsigma cut on TPC only"};
    Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 99., "Pion nsigma cut on TPC only"};
    Configurable<double> nSigmaTofPrMin{"nSigmaTofPrMin", -99, "Proton nsigma cut on TPC only"};
    Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 99., "Proton nsigma cut on TPC only"};
    Configurable<double> nSigmaTofElCombinedMax{"nSigmaTofElCombinedMax", 4., "Electron Nsigma cut on TOF combined with TPC"};
    Configurable<double> nSigmaTofPiCombinedMin{"nSigmaTofPiCombinedMin", 2.5, "Pion Nsigma cut on TOF combined with TPC"};
    Configurable<double> nSigmaTofPiCombinedMax{"nSigmaTofPiCombinedMax", 99., "Pion Nsigma cut on TOF combined with TPC"};
    Configurable<double> nSigmaTofPrCombinedMin{"nSigmaTofPrCombinedMin", -99., "Proton Nsigma cut on TOF combined with TPC"};
    Configurable<double> nSigmaTofPrCombinedMax{"nSigmaTofPrCombinedMax", 99, "Proton Nsigma cut on TOF combined with TPC"};
    // AND logic for TOF+TPC PID (as in Run2)
    Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", true, "Use AND logic for TPC and TOF PID"};
  } selectionsPid;
  Configurable<double> ptJpsiMin{"ptJpsiMin", 0., "Lower bound of J/Psi pT"};
  Configurable<double> ptJpsiMax{"ptJpsiMax", 50., "Upper bound of J/Psi pT"};
  Configurable<bool> useTrackIsGlobalTrackWoDCA{"useTrackIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for the bachelor tracks"};
  Configurable<double> ptTrackMin{"ptTrackMin", 0.5, "minimum bachelor track pT threshold (GeV/c)"};
  Configurable<double> absEtaTrackMax{"absEtaTrackMax", 0.8, "maximum bachelor track absolute eta threshold"};
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for bachelor track DCA XY pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsTrackDCA{"cutsTrackDCA", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for bachelor track"};
  // topological/kinematic cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_jpsi_to_mu_mu::vecBinsPt}, "J/Psi pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_jpsi_to_mu_mu::Cuts[0], hf_cuts_jpsi_to_mu_mu::NBinsPt, hf_cuts_jpsi_to_mu_mu::NCutVars, hf_cuts_jpsi_to_mu_mu::labelsPt, hf_cuts_jpsi_to_mu_mu::labelsCutVar}, "J/Psi candidate selection per pT bin"};
  Configurable<double> invMassWindowJpsiHad{"invMassWindowJpsiHad", 0.3, "invariant-mass window for Jpsi-Had pair preselections (GeV/c2)"};
  Configurable<double> deltaMPhiMax{"deltaMPhiMax", 0.02, "invariant-mass window for phi preselections (GeV/c2) (only for Bs->J/PsiPhi)"};
  Configurable<double> cpaMin{"cpaMin", 0., "Minimum cosine of pointing angle for B candidates"};
  Configurable<double> decLenMin{"decLenMin", 0., "Minimum decay length for B candidates"};

  // magnetic field setting from CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelper hfHelper;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;
  TrackSelectorEl selectorElectron;

  // CCDB service
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl>;
  using TracksPidWithSel = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;
  using TracksPidWithSelAndMc = soa::Join<TracksPidWithSel, aod::McTrackLabels>;
  using CollisionsWCMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using BCsInfo = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  Preslice<aod::HfCand2ProngWPid> candsJpsiPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CollisionsWCMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};
  double bz{0.};
  double invMass2JpsiHadMin{}, invMass2JpsiHadMax{};
  bool isHfCandBhadConfigFilled = false;

  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::hf_evsel::HfEventSelectionMc hfEvSelMc;

  o2::vertexing::DCAFitterN<2> df2;
  o2::vertexing::DCAFitterN<3> df3;
  o2::vertexing::DCAFitterN<4> df4;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    selectorPion.setRangePtTpc(selectionsPid.ptPidTpcMin, selectionsPid.ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-selectionsPid.nSigmaTpcPiMax, selectionsPid.nSigmaTpcPiMax);
    selectorPion.setRangeNSigmaTpcCondTof(-selectionsPid.nSigmaTpcPiCombinedMax, selectionsPid.nSigmaTpcPiCombinedMax);
    selectorPion.setRangePtTof(selectionsPid.ptPidTofMin, selectionsPid.ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-selectionsPid.nSigmaTofPiMax, selectionsPid.nSigmaTofPiMax);
    selectorPion.setRangeNSigmaTofCondTpc(-selectionsPid.nSigmaTofPiCombinedMax, selectionsPid.nSigmaTofPiCombinedMax);

    selectorProton.setRangePtTpc(selectionsPid.ptPidTpcMin, selectionsPid.ptPidTpcMax);
    selectorProton.setRangeNSigmaTpc(-selectionsPid.nSigmaTpcPrMax, selectionsPid.nSigmaTpcPrMax);
    selectorProton.setRangeNSigmaTpcCondTof(-selectionsPid.nSigmaTpcPrCombinedMax, selectionsPid.nSigmaTpcPrCombinedMax);
    selectorProton.setRangePtTof(selectionsPid.ptPidTofMin, selectionsPid.ptPidTofMax);
    selectorProton.setRangeNSigmaTof(-selectionsPid.nSigmaTofPrMax, selectionsPid.nSigmaTofPrMax);
    selectorProton.setRangeNSigmaTofCondTpc(-selectionsPid.nSigmaTofPrCombinedMax, selectionsPid.nSigmaTofPrCombinedMax);

    selectorElectron.setRangePtTpc(selectionsPid.ptPidTpcMin, selectionsPid.ptPidTpcMax);
    selectorElectron.setRangeNSigmaTpc(-selectionsPid.nSigmaTpcElMax, selectionsPid.nSigmaTpcElMax);
    selectorElectron.setRangeNSigmaTpcCondTof(-selectionsPid.nSigmaTofElCombinedMax, selectionsPid.nSigmaTofElCombinedMax);
    selectorElectron.setRangePtTof(selectionsPid.ptPidTofMin, selectionsPid.ptPidTofMax);
    selectorElectron.setRangeNSigmaTof(-selectionsPid.nSigmaTofElMax, selectionsPid.nSigmaTofElMax);
    selectorElectron.setRangeNSigmaTofCondTpc(-selectionsPid.nSigmaTofElCombinedMax, selectionsPid.nSigmaTofElCombinedMax);

    std::array<int, 4> doProcess = {doprocessJpsiKData, doprocessJpsiKMc, doprocessJpsiPhiData, doprocessJpsiPhiMc};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function can be enabled at a time, please fix your configuration!");
    }

    // Set up the histogram registry
    constexpr int kNBinsSelections = 2 + aod::SelectionStep::RecoPID;
    std::string labels[kNBinsSelections];
    labels[0] = "No selection";
    labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
    labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
    labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
    static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
    registry.add("hSelectionsJpsi", "J/Psi selection;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
      registry.get<TH2>(HIST("hSelectionsJpsi"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    constexpr int kNBinsEvents = NEvent;
    std::string labelsEvents[kNBinsEvents];
    labelsEvents[Event::Processed] = "processed";
    labelsEvents[Event::NoCharmHadPiSelected] = "without CharmHad-Pi pairs";
    labelsEvents[Event::CharmHadPiSelected] = "with CharmHad-Pi pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labelsEvents[iBin].data());
    }

    registry.add("hMassJpsi", "J/Psi mass;#it{M}_{#mu#mu} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{600, 2.8, 3.4, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtJpsi", "J/Psi #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCpaJpsi", "J/Psi cos#theta_{p};J/Psi cos#theta_{p};Counts", {HistType::kTH1F, {{200, -1., 1, "J/Psi cos#theta_{p}"}}});
    std::shared_ptr<TH1> hFitCandidatesJpsi = registry.add<TH1>("hFitCandidatesJpsi", "Jpsi candidate counter", {HistType::kTH1D, {axisCands}});
    std::shared_ptr<TH1> hFitCandidatesBPlus = registry.add<TH1>("hFitCandidatesBPlus", "hFitCandidatesBPlus candidate counter", {HistType::kTH1D, {axisCands}});
    std::shared_ptr<TH1> hFitCandidatesBS = registry.add<TH1>("hFitCandidatesBS", "hFitCandidatesBS candidate counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hFitCandidatesJpsi);
    setLabelHistoCands(hFitCandidatesBPlus);
    setLabelHistoCands(hFitCandidatesBS);
    if (doprocessJpsiKData || doprocessJpsiKMc) {
      registry.add("hPtKaon", "Kaon #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hMassJpsiKaon", "J/Psi Kaon mass;#it{M}_{J/#PsiK} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{800, 4.9, 5.7}}});
    } else if (doprocessJpsiPhiData || doprocessJpsiPhiMc) {
      registry.add("hPtPhi", "Phi #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hMassPhi", "Phi mass;#it{M}_{KK} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{200, 0.9, 1.2}}});
      registry.add("hMassJpsiPhi", "J/Psi Phi mass;#it{M}_{J/#Psi#phi} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{800, 4.9, 5.7}}});
      std::shared_ptr<TH1> hFitCandidatesPhi = registry.add<TH1>("hFitCandidatesPhi", "Phi candidate counter", {HistType::kTH1D, {axisCands}});
      setLabelHistoCands(hFitCandidatesPhi);
    }

    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
    df2.setMatCorrType(noMatCorr);

    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);
    df3.setMatCorrType(noMatCorr);

    df4.setPropagateToPCA(propagateToPCA);
    df4.setMaxR(maxR);
    df4.setMaxDZIni(maxDZIni);
    df4.setMinParamChange(minParamChange);
    df4.setMinRelChi2Change(minRelChi2Change);
    df4.setUseAbsDCA(useAbsDCA);
    df4.setWeightedFinalPCA(useWeightedFinalPCA);
    df4.setMatCorrType(noMatCorr);

    // Configure CCDB access
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    if (doprocessJpsiKData || doprocessJpsiKMc) {
      invMass2JpsiHadMin = (MassBPlus - invMassWindowJpsiHad) * (MassBPlus - invMassWindowJpsiHad);
      invMass2JpsiHadMax = (MassBPlus + invMassWindowJpsiHad) * (MassBPlus + invMassWindowJpsiHad);
    } else if (doprocessJpsiPhiData || doprocessJpsiPhiMc) {
      invMass2JpsiHadMin = (MassBS - invMassWindowJpsiHad) * (MassBS - invMassWindowJpsiHad);
      invMass2JpsiHadMax = (MassBS + invMassWindowJpsiHad) * (MassBS + invMassWindowJpsiHad);
    }

    // init HF event selection helper
    hfEvSel.init(registry);
    if (doprocessJpsiKMc || doprocessJpsiPhiMc) {
      const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
      for (const DeviceSpec& device : workflows.devices) {
        if (device.name == "hf-data-creator-jpsi-had-reduced") {
          // init HF event selection helper
          hfEvSelMc.init(device, registry);
          break;
        }
      }
    }
  }

  /// Topological cuts
  /// \param candidate is candidate
  /// \param trackPos is the positive track
  /// \param trackNeg is the negative track
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionTopol(const T1& candidate, const T2& trackPos, const T2& trackNeg)
  {
    auto candpT = candidate.pt();
    auto candInvMass = runJpsiToee ? hfHelper.invMassJpsiToEE(candidate) : hfHelper.invMassJpsiToMuMu(candidate);
    auto pseudoPropDecLen = candidate.decayLengthXY() * candInvMass / candpT;
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptJpsiMin || candpT >= ptJpsiMax) {
      return false;
    }

    // cut on μ+ μ− (e+e−) invariant mass
    if (std::abs(candInvMass - o2::constants::physics::MassJPsi) > cuts->get(pTBin, "m")) {
      return false;
    }

    // cut on daughter pT (same cut used for both channels)
    if (trackNeg.pt() < cuts->get(pTBin, "pT mu") || trackPos.pt() < cuts->get(pTBin, "pT mu")) {
      return false;
    }

    // decay length
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }

    // decay length in XY plane
    if (candidate.decayLengthXY() < cuts->get(pTBin, "decay length xy")) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() < cuts->get(pTBin, "cpa")) {
      return false;
    }

    // cosine of pointing angle XY
    if (candidate.cpaXY() < cuts->get(pTBin, "cpa xy")) {
      return false;
    }

    // product of daughter impact parameters
    if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0xd0")) {
      return false;
    }

    // pseudoproper decay length
    if (pseudoPropDecLen < cuts->get(pTBin, "pseudoprop. decay length")) {
      return false;
    }

    return true;
  }

  /// Kaon selection (J/Psi K+ <-- B+)
  /// \param track is the considered track
  /// \param trackParCov is the track parametrisation
  /// \param dca is the 2-D array with track DCAs
  /// \param jPsiDautracks J/Psi daughter tracks
  /// \return true if track passes all cuts
  template <typename T1, typename T2, typename T3>
  bool isTrackSelected(const T1& track, const T2& trackParCov, const T3& dca, const std::vector<T1>& jPsiDautracks)
  {
    // check isGlobalTrackWoDCA status for kaons if wanted
    if (useTrackIsGlobalTrackWoDCA && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    // minimum pT, eta, and DCA selection
    if (trackParCov.getPt() < ptTrackMin || std::abs(trackParCov.getEta()) > absEtaTrackMax || !isSelectedTrackDCA(trackParCov, dca, binsPtTrack, cutsTrackDCA)) {
      return false;
    }
    // reject kaons that are J/Psi daughters
    for (const auto& trackJpsi : jPsiDautracks) {
      if (track.globalIndex() == trackJpsi.globalIndex()) {
        return false;
      }
    }

    return true;
  }

  template <typename T1>
  bool isSelectedJpsiDauPid(const T1& track)
  {
    int pidPion = -1;
    int pidProton = -1;
    int pidElectron = -1;

    if (selectionsPid.usePidTpcAndTof) {
      pidPion = selectorPion.statusTpcAndTof(track, track.tpcNSigmaPi(), track.tofNSigmaPi());
      pidProton = selectorProton.statusTpcAndTof(track, track.tpcNSigmaPr(), track.tofNSigmaPr());
      pidElectron = selectorElectron.statusTpcAndTof(track, track.tpcNSigmaEl(), track.tofNSigmaEl());
    } else {
      pidPion = selectorPion.statusTpcOrTof(track, track.tpcNSigmaPi(), track.tofNSigmaPi());
      pidProton = selectorProton.statusTpcOrTof(track, track.tpcNSigmaPr(), track.tofNSigmaPr());
      pidElectron = selectorElectron.statusTpcOrTof(track, track.tpcNSigmaEl(), track.tofNSigmaEl());
    }

    if (pidPion == TrackSelectorPID::Rejected ||
        pidProton == TrackSelectorPID::Rejected ||
        pidElectron == TrackSelectorPID::Rejected) {
      return false;
    }
    return true;
  }

  /// B meson preselections
  /// \param momentum is the B meson momentum
  /// \param secondaryVertex is the reconstructed secondary vertex
  /// \param collision is the reconstructed collision
  template <typename T1, typename T2, typename T3>
  bool isBSelected(const T1& momentum, const T2& secondaryVertex, const T3& collision)
  {
    // B candidate CPA
    if (RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, momentum) < cpaMin) {
      return false;
    }

    // B candidate decay length
    if (RecoDecay::distance(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex) < decLenMin) {
      return false;
    }

    return true;
  }

  /// Checks if the B meson is associated with a different collision than the one it was generated in
  /// \param particleMother is the mother particle
  /// \param collision is the reconstructed collision
  /// \param indexCollisionMaxNumContrib is the index of the collision associated with a given MC collision with the largest number of contributors.
  /// \param flagWrongCollision is the flag indicating if whether the associated collision is incorrect.
  template <typename PParticle, typename CColl>
  void checkWrongCollision(const PParticle& particleMother,
                           const CColl& collision,
                           const int64_t& indexCollisionMaxNumContrib,
                           int8_t& flagWrongCollision)
  {

    if (particleMother.mcCollision().globalIndex() != collision.mcCollisionId()) {
      flagWrongCollision = WrongCollisionType::WrongAssociation;
    } else {
      if (collision.globalIndex() != indexCollisionMaxNumContrib) {
        flagWrongCollision = WrongCollisionType::SplitCollision;
      }
    }
  }

  /// Function for filling MC reco information in the tables
  /// \param particlesMc is the table with MC particles
  /// \param vecDaughtersB is the vector with all daughter tracks (Jpsi daughters in first position)
  /// \param indexHfCandJpsi is the index of the Jpsi candidate
  /// \param selectedTracksBach is the map with the indices of selected bachelor pion tracks
  template <uint8_t DecChannel, typename CColl, typename PParticles, typename TTrack>
  void fillMcRecoInfo(CColl const& collision,
                      PParticles const& particlesMc,
                      std::vector<TTrack> const& vecDaughtersB,
                      const int64_t indexHfCandJpsi,
                      std::array<std::map<int64_t, int64_t>, 2> selectedTracksBach,
                      const int64_t indexCollisionMaxNumContrib)
  {

    // we check the MC matching to be stored
    int8_t sign{0}, flag{0}, channel{0};
    int8_t flagWrongCollision{WrongCollisionType::None};
    int8_t debug{0};
    float motherPt{-1.f};

    if constexpr (DecChannel == DecayChannel::BplusToJpsiK) {
      // B+ → J/Psi K+ → (µ+µ-) K+
      int indexRec = -1;
      if (!runJpsiToee) {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{-kMuonMinus, +kMuonMinus, +kKPlus}, true, &sign, 3);
      } else {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2]}, Pdg::kBPlus, std::array{-kElectron, +kElectron, +kKPlus}, true, &sign, 3);
      }
      if (indexRec > -1) {
        // J/Psi → µ+µ-
        if (!runJpsiToee) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kMuonMinus, +kMuonMinus}, true);
        } else {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kElectron, +kElectron}, true);
        }
        if (indexRec > -1) {
          flag = sign * o2::hf_decay::hf_cand_beauty::BplusToJpsiK;
        } else {
          debug = 1;
          LOGF(debug, "B+ decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBPlus, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }
      rowHfJpsiKMcRecReduced(indexHfCandJpsi, selectedTracksBach[0][vecDaughtersB.back().globalIndex()], flag, channel, flagWrongCollision, debug, motherPt);
    } else if constexpr (DecChannel == DecayChannel::BsToJpsiPhi) {
      // Bs → J/Psi phi → (µ+µ-) (K+K-)
      int indexRec = -1;
      if (!runJpsiToee) {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kMuonMinus, +kMuonMinus, +kKPlus, -kKPlus}, true, &sign, 4);
      } else {
        indexRec = RecoDecay::getMatchedMCRec<true, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1], vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kBS, std::array{-kElectron, +kElectron, +kKPlus, -kKPlus}, true, &sign, 4);
      }
      if (indexRec > -1) {
        // J/Psi → µ+µ-
        if (!runJpsiToee) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kMuonMinus, +kMuonMinus}, true, &sign, 1);
        } else {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[0], vecDaughtersB[1]}, Pdg::kJPsi, std::array{-kElectron, +kElectron}, true, &sign, 1);
        }
        if (indexRec > -1) {
          flag = sign * o2::hf_decay::hf_cand_beauty::BsToJpsiKK;
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particlesMc, std::array{vecDaughtersB[2], vecDaughtersB[3]}, Pdg::kPhi, std::array{-kKPlus, +kKPlus}, true, &sign, 1);
          if (indexRec > -1) {
            channel = o2::hf_decay::hf_cand_beauty::BsToJpsiPhi;
          } else {
            debug = 1;
            LOGF(debug, "Bs decays in the expected final state but the condition on the phi intermediate state is not fulfilled");
          }
        } else {
          debug = 1;
          LOGF(debug, "Bs decays in the expected final state but the condition on the J/Psi intermediate state is not fulfilled");
        }

        auto indexMother = RecoDecay::getMother(particlesMc, vecDaughtersB.back().template mcParticle_as<PParticles>(), Pdg::kBS, true);
        if (indexMother >= 0) {
          auto particleMother = particlesMc.rawIteratorAt(indexMother);
          motherPt = particleMother.pt();
          checkWrongCollision(particleMother, collision, indexCollisionMaxNumContrib, flagWrongCollision);
        }
      }
      rowHfJpsiPhiMcRecReduced(indexHfCandJpsi, selectedTracksBach[0][vecDaughtersB.back().globalIndex()], selectedTracksBach[1][vecDaughtersB.back().globalIndex()], flag, channel, flagWrongCollision, debug, motherPt);
    }
  }

  /// Calculates the index of the collision with the maximum number of contributions.
  ///\param collisions are the collisions to search through.
  ///\return The index of the collision with the maximum number of contributions.
  template <typename CColl>
  int64_t getIndexCollisionMaxNumContrib(const CColl& collisions)
  {
    unsigned maxNumContrib = 0;
    int64_t indexCollisionMaxNumContrib = -1;
    for (const auto& collision : collisions) {
      if (collision.numContrib() > maxNumContrib) {
        maxNumContrib = collision.numContrib();
        indexCollisionMaxNumContrib = collision.globalIndex();
      }
    }
    return indexCollisionMaxNumContrib;
  }

  template <uint8_t DecChannel>
  void runMcGen(aod::McCollision const& mcCollision,
                aod::McParticles const& particlesMc,
                CollisionsWCMcLabels const& collisions,
                BCsInfo const&)
  {
    // Check event selection
    float centDummy{-1.f}, centFT0C{-1.f}, centFT0M{-1.f};
    const auto collSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
    auto hfRejMap = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, o2::hf_centrality::CentralityEstimator::None>(mcCollision, collSlice, centDummy);
    if (skipRejectedCollisions && hfRejMap != 0) {
      return;
    }

    const auto mcParticlesPerMcColl = particlesMc.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());

    // Match generated particles.
    for (const auto& particle : mcParticlesPerMcColl) {
      int8_t sign{0}, flag{0}, channel{0};
      if constexpr (DecChannel == DecayChannel::BplusToJpsiK) {
        // B+ → J/Psi K+ → (µ+µ-) K+
        if (RecoDecay::isMatchedMCGen<false>(particlesMc, particle, Pdg::kBPlus, std::array{static_cast<int>(Pdg::kJPsi), +kKPlus}, true, &sign)) {
          // Match J/Psi -> µ+µ-
          auto candJpsiMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking J/Psi -> µ+µ-");
          if (!runJpsiToee) {
            if (RecoDecay::isMatchedMCGen(particlesMc, candJpsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kMuonMinus, +kMuonMinus}, true)) {
              flag = sign * o2::hf_decay::hf_cand_beauty::BplusToJpsiK;
            }
          } else { // debug
            // Printf("Checking J/Psi -> e+e-");
            if (RecoDecay::isMatchedMCGen(particlesMc, candJpsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kElectron, +kElectron}, true)) {
              flag = sign * o2::hf_decay::hf_cand_beauty::BplusToJpsiK;
            }
          }
        }

        // save information for B+ task
        if (std::abs(flag) != o2::hf_decay::hf_cand_beauty::BplusToJpsiK) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassBPlus);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfBpMcGenReduced(flag, channel, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      } else if constexpr (DecChannel == DecayChannel::BsToJpsiPhi) {
        // Bs → J/Psi phi → (µ+µ-) (K+K-)
        if (RecoDecay::isMatchedMCGen<true>(particlesMc, particle, Pdg::kBS, std::array{static_cast<int>(Pdg::kJPsi), +kKPlus, -kKPlus}, true, &sign, 2)) {
          // Match J/Psi -> µ+µ- and phi -> K+K-
          auto candJpsiMC = particlesMc.rawIteratorAt(particle.daughtersIds().front());
          auto candPhiMC = particlesMc.rawIteratorAt(particle.daughtersIds().back());
          // Printf("Checking J/Psi -> µ+µ- and phi -> K+K-");
          if (runJpsiToee && RecoDecay::isMatchedMCGen(particlesMc, candJpsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kElectron, +kElectron}, true)) {
            flag = sign * o2::hf_decay::hf_cand_beauty::BsToJpsiKK;
          } else if (!runJpsiToee && RecoDecay::isMatchedMCGen(particlesMc, candJpsiMC, static_cast<int>(Pdg::kJPsi), std::array{-kMuonMinus, +kMuonMinus}, true)) {
            flag = sign * o2::hf_decay::hf_cand_beauty::BsToJpsiKK;
          }
          // Check phi -> K+K-
          if (RecoDecay::isMatchedMCGen(particlesMc, candPhiMC, static_cast<int>(Pdg::kPhi), std::array{-kKPlus, +kKPlus}, true)) {
            channel = o2::hf_decay::hf_cand_beauty::BsToJpsiPhi;
          }
        }

        // save information for Bs task
        if (std::abs(flag) != o2::hf_decay::hf_cand_beauty::BsToJpsiKK) {
          continue;
        }

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), MassBPlus);
        auto etaParticle = particle.eta();

        std::array<float, 2> ptProngs{};
        std::array<float, 2> yProngs{};
        std::array<float, 2> etaProngs{};
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }
        rowHfBsMcGenReduced(flag, channel, ptParticle, yParticle, etaParticle,
                            ptProngs[0], yProngs[0], etaProngs[0],
                            ptProngs[1], yProngs[1], etaProngs[1], hfRejMap, centFT0C, centFT0M);
      }
    } // gen
  }

  // Jpsi candidate selection
  template <bool DoMc, uint8_t DecChannel, typename Coll, typename JpsiCands, typename TTracks, typename PParticles, typename BBCs>
  void runDataCreation(Coll const& collision,
                       JpsiCands const& candsJpsi,
                       aod::TrackAssoc const& trackIndices,
                       TTracks const&,
                       PParticles const& particlesMc,
                       uint64_t const& indexCollisionMaxNumContrib,
                       BBCs const&)
  {

    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    float centrality = -1.f;
    const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    if (skipRejectedCollisions && hfRejMap != 0) {
      return;
    }

    // helpers for ReducedTables filling
    int const indexHfReducedCollision = hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the track.globalIndex() and
    // the value is the track index in the table of the selected tracks
    std::map<int64_t, int64_t> selectedTracksBach;
    std::map<int64_t, int64_t> selectedTracksBach2; // for the second daughter (for B0 and Bs)

    bool fillHfReducedCollision = false;

    auto primaryVertex = getPrimaryVertex(collision);

    // Set the magnetic field from ccdb.
    // The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
    // but this is not true when running on Run2 data/MC already converted into AO2Ds.
    auto bc = collision.template bc_as<BBCs>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    df2.setBz(bz);
    df3.setBz(bz);
    df4.setBz(bz);

    auto thisCollId = collision.globalIndex();
    // looping over 2-prong candidates
    for (const auto& candidate : candsJpsi) {

      // Apply the selections on the J/Psi candidates
      registry.fill(HIST("hSelectionsJpsi"), 1, candidate.pt());

      if (!(candidate.hfflag() & (1 << (runJpsiToee ? aod::hf_cand_2prong::DecayType::JpsiToEE : aod::hf_cand_2prong::DecayType::JpsiToMuMu)))) {
        continue;
      }
      registry.fill(HIST("hSelectionsJpsi"), 2 + aod::SelectionStep::RecoSkims, candidate.pt());

      auto trackPos = candidate.template prong0_as<TTracks>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<TTracks>(); // negative daughter

      auto trackPosParCov = getTrackParCov(trackPos);
      auto trackNegParCov = getTrackParCov(trackNeg);

      std::vector<typename TTracks::iterator> jPsiDauTracks{trackPos, trackNeg};

      auto dca0 = o2::dataformats::DCA(jPsiDauTracks[0].dcaXY(), jPsiDauTracks[0].dcaZ(), jPsiDauTracks[0].cYY(), jPsiDauTracks[0].cZY(), jPsiDauTracks[0].cZZ());
      auto dca1 = o2::dataformats::DCA(jPsiDauTracks[1].dcaXY(), jPsiDauTracks[1].dcaZ(), jPsiDauTracks[1].cYY(), jPsiDauTracks[1].cZY(), jPsiDauTracks[1].cZZ());

      // repropagate tracks to this collision if needed
      if (jPsiDauTracks[0].collisionId() != thisCollId) {
        trackPosParCov.propagateToDCA(primaryVertex, bz, &dca0);
      }

      if (jPsiDauTracks[1].collisionId() != thisCollId) {
        trackNegParCov.propagateToDCA(primaryVertex, bz, &dca1);
      }

      // ---------------------------------
      // reconstruct J/Psi candidate secondary vertex
      o2::track::TrackParCov const trackParCovJpsi{}; // FIXME: unused
      std::array<float, 3> const pVecJpsi{};          // FIXME: unused
      registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::BeforeFit);
      try {
        if (df2.process(trackPosParCov, trackNegParCov) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::Fail);
        continue;
      }
      registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::FitOk);

      // topological selection
      if (!selectionTopol(candidate, trackPos, trackNeg)) {
        continue;
      }
      registry.fill(HIST("hSelectionsJpsi"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());

      // PID selection
      if (!isSelectedJpsiDauPid(trackPos) || !isSelectedJpsiDauPid(trackNeg)) {
        continue;
      }
      registry.fill(HIST("hSelectionsJpsi"), 2 + aod::SelectionStep::RecoPID, candidate.pt());

      int const indexHfCandJpsi = hfJpsi.lastIndex() + 1;
      float const invMassJpsi = runJpsiToee ? hfHelper.invMassJpsiToEE(candidate) : hfHelper.invMassJpsiToMuMu(candidate);
      registry.fill(HIST("hMassJpsi"), invMassJpsi);
      registry.fill(HIST("hPtJpsi"), candidate.pt());
      registry.fill(HIST("hCpaJpsi"), candidate.cpa());

      bool fillHfCandJpsi = false;

      // TODO: add single track information (min eta, min ITS/TPC clusters, etc.)
      double invMass2JpsiHad{0.};
      for (const auto& trackId : trackIndices) {
        auto trackBach = trackId.template track_as<TTracks>();

        // apply selections on bachelor tracks
        auto trackParCovBach = getTrackParCov(trackBach);
        std::array<float, 2> dcaBach{trackBach.dcaXY(), trackBach.dcaZ()};
        std::array<float, 3> pVecBach = trackBach.pVector();
        if (trackBach.collisionId() != thisCollId) {
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovBach, 2.f, noMatCorr, &dcaBach);
          getPxPyPz(trackParCovBach, pVecBach);
        }

        // apply selections on bachelor tracks
        if (!isTrackSelected(trackBach, trackParCovBach, dcaBach, jPsiDauTracks)) {
          continue;
        }

        if constexpr (DecChannel == DecayChannel::BplusToJpsiK) {
          registry.fill(HIST("hPtKaon"), trackParCovBach.getPt());
          // compute invariant mass square and apply selection
          invMass2JpsiHad = RecoDecay::m2(std::array{pVecJpsi, pVecBach}, std::array{MassJPsi, MassKPlus});
          if ((invMass2JpsiHad < invMass2JpsiHadMin) || (invMass2JpsiHad > invMass2JpsiHadMax)) {
            continue;
          }
          registry.fill(HIST("hMassJpsiKaon"), std::sqrt(invMass2JpsiHad));

          registry.fill(HIST("hFitCandidatesBPlus"), SVFitting::BeforeFit);
          try {
            if (df3.process(trackPosParCov, trackNegParCov, trackParCovBach) == 0) {
              continue;
            }
          } catch (const std::runtime_error& error) {
            LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
            registry.fill(HIST("hFitCandidatesBPlus"), SVFitting::Fail);
            continue;
          }
          registry.fill(HIST("hFitCandidatesBPlus"), SVFitting::FitOk);

          o2::track::TrackParCov trackParCovBPlus{};
          std::array<float, 3> pVecBPlus{}, pVec0{}, pVec1{}, pVec2{};

          auto secondaryVertexBPlus = df3.getPCACandidate();
          df3.propagateTracksToVertex();
          df3.getTrack(0).getPxPyPzGlo(pVec0);
          df3.getTrack(1).getPxPyPzGlo(pVec1);
          df3.getTrack(2).getPxPyPzGlo(pVec2);
          pVecBPlus = RecoDecay::pVec(pVec0, pVec1, pVec2);
          trackParCovBPlus = df3.createParentTrackParCov();
          trackParCovBPlus.setAbsCharge(0); // to be sure

          if (!isBSelected(pVecBPlus, secondaryVertexBPlus, collision)) {
            continue;
          }

          // fill Kaon tracks table
          // if information on track already stored, go to next track
          if (!selectedTracksBach.count(trackBach.globalIndex())) {
            hfTrackLfDau0(trackBach.globalIndex(), indexHfReducedCollision,
                          trackParCovBach.getX(), trackParCovBach.getAlpha(),
                          trackParCovBach.getY(), trackParCovBach.getZ(), trackParCovBach.getSnp(),
                          trackParCovBach.getTgl(), trackParCovBach.getQ2Pt(),
                          trackBach.itsNCls(), trackBach.tpcNClsCrossedRows(), trackBach.tpcChi2NCl(), trackBach.itsChi2NCl(),
                          trackBach.hasTPC(), trackBach.hasTOF(),
                          trackBach.tpcNSigmaPi(), trackBach.tofNSigmaPi(),
                          trackBach.tpcNSigmaKa(), trackBach.tofNSigmaKa(),
                          trackBach.tpcNSigmaPr(), trackBach.tofNSigmaPr());
            hfTrackCovLfDau0(trackParCovBach.getSigmaY2(), trackParCovBach.getSigmaZY(), trackParCovBach.getSigmaZ2(),
                             trackParCovBach.getSigmaSnpY(), trackParCovBach.getSigmaSnpZ(),
                             trackParCovBach.getSigmaSnp2(), trackParCovBach.getSigmaTglY(), trackParCovBach.getSigmaTglZ(),
                             trackParCovBach.getSigmaTglSnp(), trackParCovBach.getSigmaTgl2(),
                             trackParCovBach.getSigma1PtY(), trackParCovBach.getSigma1PtZ(), trackParCovBach.getSigma1PtSnp(),
                             trackParCovBach.getSigma1PtTgl(), trackParCovBach.getSigma1Pt2());
            // add trackBach.globalIndex() to a list
            // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another Jpsi candidate
            // and keep track of their index in hfTrackLfDau0 for McRec purposes
            selectedTracksBach[trackBach.globalIndex()] = hfTrackLfDau0.lastIndex();
          }

          if constexpr (DoMc) {
            std::vector<typename TTracks::iterator> beautyHadDauTracks{};
            beautyHadDauTracks.reserve(jPsiDauTracks.size());
            for (const auto& track : jPsiDauTracks) {
              beautyHadDauTracks.push_back(track);
            }
            beautyHadDauTracks.push_back(trackBach);
            fillMcRecoInfo<DecayChannel::BplusToJpsiK>(collision, particlesMc, beautyHadDauTracks, indexHfCandJpsi, std::array<std::map<int64_t, int64_t>, 2>{selectedTracksBach}, indexCollisionMaxNumContrib);
          }
          fillHfCandJpsi = true;
        } else if constexpr (DecChannel == DecayChannel::BsToJpsiPhi) {
          for (auto trackBachId2 = trackId + 1; trackBachId2 != trackIndices.end(); ++trackBachId2) {
            auto trackBach2 = trackBachId2.template track_as<TTracks>();
            auto trackBach2ParCov = getTrackParCov(trackBach2);

            std::array<float, 2> dcaBach2{trackBach2.dcaXY(), trackBach2.dcaZ()};
            std::array<float, 3> pVecBach2 = trackBach2.pVector();
            if (trackBach2.collisionId() != thisCollId) {
              o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackBach2ParCov, 2.f, noMatCorr, &dcaBach2);
              getPxPyPz(trackBach2ParCov, pVecBach2);
            }

            // apply selections on bachelor tracks
            if (!isTrackSelected(trackBach2, trackBach2ParCov, dcaBach2, jPsiDauTracks)) {
              continue;
            }
            std::array<float, 3> pVec2{trackBach.pVector()}, pVec3{trackBach2.pVector()};
            auto invMassPhi = RecoDecay::m(std::array{pVec2, pVec3}, std::array{MassKPlus, MassKPlus});

            if (std::abs(invMassPhi - MassPhi) > deltaMPhiMax) {
              continue;
            }

            // ---------------------------------
            // reconstruct Bs candidate secondary vertex

            registry.fill(HIST("hFitCandidatesBS"), SVFitting::BeforeFit);
            try {
              if (df4.process(trackPosParCov, trackNegParCov, trackParCovBach, trackBach2ParCov) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
              registry.fill(HIST("hFitCandidatesBS"), SVFitting::Fail);
              continue;
            }
            registry.fill(HIST("hFitCandidatesBS"), SVFitting::FitOk);

            o2::track::TrackParCov trackParCovBS{};
            std::array<float, 3> pVecBS{}, pVec0{}, pVec1{}, pVecPhi{};

            auto secondaryVertexBS = df4.getPCACandidate();
            df4.propagateTracksToVertex();
            df4.getTrack(0).getPxPyPzGlo(pVec0);
            df4.getTrack(1).getPxPyPzGlo(pVec1);
            df4.getTrack(2).getPxPyPzGlo(pVec2);
            df4.getTrack(3).getPxPyPzGlo(pVec3);
            pVecBS = RecoDecay::pVec(pVec0, pVec1, pVec2, pVec3);
            pVecPhi = RecoDecay::pVec(pVec2, pVec3);
            trackParCovBS = df4.createParentTrackParCov();
            trackParCovBS.setAbsCharge(0); // to be sure

            if (!isBSelected(pVecBS, secondaryVertexBS, collision)) {
              continue;
            }

            registry.fill(HIST("hPtPhi"), RecoDecay::pt(pVecBach, pVecBach2));
            registry.fill(HIST("hMassPhi"), RecoDecay::m(std::array{pVecBach, pVecBach2}, std::array{MassKPlus, MassKPlus}));
            invMass2JpsiHad = RecoDecay::m2(std::array{pVecJpsi, pVecPhi}, std::array{MassJPsi, MassPhi});
            if ((invMass2JpsiHad < invMass2JpsiHadMin) || (invMass2JpsiHad > invMass2JpsiHadMax)) {
              continue;
            }
            registry.fill(HIST("hMassJpsiPhi"), std::sqrt(invMass2JpsiHad));

            // fill daughter tracks table
            // if information on track already stored, go to next track
            if (!selectedTracksBach.count(trackBach.globalIndex())) {
              hfTrackLfDau0(trackBach.globalIndex(), indexHfReducedCollision,
                            trackParCovBach.getX(), trackParCovBach.getAlpha(),
                            trackParCovBach.getY(), trackParCovBach.getZ(), trackParCovBach.getSnp(),
                            trackParCovBach.getTgl(), trackParCovBach.getQ2Pt(),
                            trackBach.itsNCls(), trackBach.tpcNClsCrossedRows(), trackBach.tpcChi2NCl(), trackBach.itsChi2NCl(),
                            trackBach.hasTPC(), trackBach.hasTOF(),
                            trackBach.tpcNSigmaPi(), trackBach.tofNSigmaPi(),
                            trackBach.tpcNSigmaKa(), trackBach.tofNSigmaKa(),
                            trackBach.tpcNSigmaPr(), trackBach.tofNSigmaPr());
              hfTrackCovLfDau0(trackParCovBach.getSigmaY2(), trackParCovBach.getSigmaZY(), trackParCovBach.getSigmaZ2(),
                               trackParCovBach.getSigmaSnpY(), trackParCovBach.getSigmaSnpZ(),
                               trackParCovBach.getSigmaSnp2(), trackParCovBach.getSigmaTglY(), trackParCovBach.getSigmaTglZ(),
                               trackParCovBach.getSigmaTglSnp(), trackParCovBach.getSigmaTgl2(),
                               trackParCovBach.getSigma1PtY(), trackParCovBach.getSigma1PtZ(), trackParCovBach.getSigma1PtSnp(),
                               trackParCovBach.getSigma1PtTgl(), trackParCovBach.getSigma1Pt2());
              // add trackBach.globalIndex() to a list
              // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another Jpsi candidate
              // and keep track of their index in hfTrackLfDau0 for McRec purposes
              selectedTracksBach[trackBach.globalIndex()] = hfTrackLfDau0.lastIndex();
            }

            // fill daughter tracks table
            // if information on track already stored, go to next track
            if (!selectedTracksBach2.count(trackBach2.globalIndex())) {
              hfTrackLfDau1(trackBach2.globalIndex(), indexHfReducedCollision,
                            trackBach2ParCov.getX(), trackBach2ParCov.getAlpha(),
                            trackBach2ParCov.getY(), trackBach2ParCov.getZ(), trackBach2ParCov.getSnp(),
                            trackBach2ParCov.getTgl(), trackBach2ParCov.getQ2Pt(),
                            trackBach2.itsNCls(), trackBach2.tpcNClsCrossedRows(), trackBach2.tpcChi2NCl(), trackBach2.itsChi2NCl(),
                            trackBach2.hasTPC(), trackBach2.hasTOF(),
                            trackBach2.tpcNSigmaPi(), trackBach2.tofNSigmaPi(),
                            trackBach2.tpcNSigmaKa(), trackBach2.tofNSigmaKa(),
                            trackBach2.tpcNSigmaPr(), trackBach2.tofNSigmaPr());
              hfTrackCovLfDau1(trackBach2ParCov.getSigmaY2(), trackBach2ParCov.getSigmaZY(), trackBach2ParCov.getSigmaZ2(),
                               trackBach2ParCov.getSigmaSnpY(), trackBach2ParCov.getSigmaSnpZ(),
                               trackBach2ParCov.getSigmaSnp2(), trackBach2ParCov.getSigmaTglY(), trackBach2ParCov.getSigmaTglZ(),
                               trackBach2ParCov.getSigmaTglSnp(), trackBach2ParCov.getSigmaTgl2(),
                               trackBach2ParCov.getSigma1PtY(), trackBach2ParCov.getSigma1PtZ(), trackBach2ParCov.getSigma1PtSnp(),
                               trackBach2ParCov.getSigma1PtTgl(), trackBach2ParCov.getSigma1Pt2());
              // add trackBach2.globalIndex() to a list
              // to keep memory of the pions filled in the table and avoid refilling them if they are paired to another Jpsi candidate
              // and keep track of their index in hfTrackLfDau1 for McRec purposes
              selectedTracksBach2[trackBach2.globalIndex()] = hfTrackLfDau1.lastIndex();
            }

            if constexpr (DoMc) {
              std::vector<typename TTracks::iterator> beautyHadDauTracks{};
              beautyHadDauTracks.reserve(jPsiDauTracks.size());
              for (const auto& track : jPsiDauTracks) {
                beautyHadDauTracks.push_back(track);
              }
              beautyHadDauTracks.push_back(trackBach);
              fillMcRecoInfo<DecayChannel::BplusToJpsiK>(collision, particlesMc, beautyHadDauTracks, indexHfCandJpsi, std::array<std::map<int64_t, int64_t>, 2>{selectedTracksBach, selectedTracksBach2}, indexCollisionMaxNumContrib);
            }
            fillHfCandJpsi = true;
          }
        }
      } // kaon loop
      if (fillHfCandJpsi) { // fill Jpsi table only once per Jpsi candidate
        double invMassJpsi{0.};
        if (runJpsiToee) {
          invMassJpsi = hfHelper.invMassJpsiToEE(candidate);
        } else {
          invMassJpsi = hfHelper.invMassJpsiToMuMu(candidate);
        }
        hfJpsi(trackPos.globalIndex(), trackNeg.globalIndex(),
               indexHfReducedCollision,
               candidate.xSecondaryVertex(), candidate.ySecondaryVertex(), candidate.zSecondaryVertex(),
               invMassJpsi,
               trackPos.itsNCls(), trackPos.tpcNClsCrossedRows(), trackPos.tpcChi2NCl(), trackPos.itsChi2NCl(),
               trackNeg.itsNCls(), trackNeg.tpcNClsCrossedRows(), trackNeg.tpcChi2NCl(), trackNeg.itsChi2NCl(),
               trackPosParCov.getX(), trackNegParCov.getX(),
               trackPosParCov.getY(), trackNegParCov.getY(),
               trackPosParCov.getZ(), trackNegParCov.getZ(),
               trackPosParCov.getAlpha(), trackNegParCov.getAlpha(),
               trackPosParCov.getSnp(), trackNegParCov.getSnp(),
               trackPosParCov.getTgl(), trackNegParCov.getTgl(),
               trackPosParCov.getQ2Pt(), trackNegParCov.getQ2Pt()); // Q/pT
        hfRedJpsiCov(trackPosParCov.getSigmaY2(), trackNegParCov.getSigmaY2(),
                     trackPosParCov.getSigmaZY(), trackNegParCov.getSigmaZY(),
                     trackPosParCov.getSigmaZ2(), trackNegParCov.getSigmaZ2(),
                     trackPosParCov.getSigmaSnpY(), trackNegParCov.getSigmaSnpY(),
                     trackPosParCov.getSigmaSnpZ(), trackNegParCov.getSigmaSnpZ(),
                     trackPosParCov.getSigmaSnp2(), trackNegParCov.getSigmaSnp2(),
                     trackPosParCov.getSigmaTglY(), trackNegParCov.getSigmaTglY(),
                     trackPosParCov.getSigmaTglZ(), trackNegParCov.getSigmaTglZ(),
                     trackPosParCov.getSigmaTglSnp(), trackNegParCov.getSigmaTglSnp(),
                     trackPosParCov.getSigmaTgl2(), trackNegParCov.getSigmaTgl2(),
                     trackPosParCov.getSigma1PtY(), trackNegParCov.getSigma1PtY(),
                     trackPosParCov.getSigma1PtZ(), trackNegParCov.getSigma1PtZ(),
                     trackPosParCov.getSigma1PtSnp(), trackNegParCov.getSigma1PtSnp(),
                     trackPosParCov.getSigma1PtTgl(), trackNegParCov.getSigma1PtTgl(),
                     trackPosParCov.getSigma1Pt2(), trackNegParCov.getSigma1Pt2());
        fillHfReducedCollision = true;
      }
    } // candsJpsi loop

    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoCharmHadPiSelected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::CharmHadPiSelected);
    // fill collision table if it contains a J/Psi K pair at minimum
    hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
    hfReducedCollExtra(collision.covXX(), collision.covXY(), collision.covYY(),
                       collision.covXZ(), collision.covYZ(), collision.covZZ());
    // hfReducedCollCentrality(collision.centFT0C(), collision.centFT0M(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange()); // TODO: add
    // if constexpr (withQvec) {
    //   hfReducedQvector(collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C(),
    //                    collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A(),
    //                    collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M(),
    //                    collision.qvecTPCposRe(), collision.qvecTPCposIm(), collision.nTrkTPCpos(),
    //                    collision.qvecTPCnegRe(), collision.qvecTPCnegIm(), collision.nTrkTPCneg(),
    //                    collision.qvecTPCallRe(), collision.qvecTPCallIm(), collision.nTrkTPCall());
    // }
  }

  void processJpsiKData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        aod::HfCand2ProngWPid const& candsJpsi,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSel const& tracks,
                        aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(invMassWindowJpsiHad.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, DecayChannel::BplusToJpsiK>(collision, candsJpsiThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorJpsiHadReduced, processJpsiKData, "Process J/Psi K without MC info", true);

  void processJpsiPhiData(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                          aod::HfCand2ProngWPid const& candsJpsi,
                          aod::TrackAssoc const& trackIndices,
                          TracksPidWithSel const& tracks,
                          aod::BCsWithTimestamps const& bcs)
  {
    // store configurables needed for B0 workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(invMassWindowJpsiHad.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, DecayChannel::BsToJpsiPhi>(collision, candsJpsiThisColl, trackIdsThisCollision, tracks, tracks, -1, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorJpsiHadReduced, processJpsiPhiData, "Process J/Psi phi without MC info", false);

  void processJpsiKMc(CollisionsWCMcLabels const& collisions,
                      aod::HfCand2ProngWPid const& candsJpsi,
                      aod::TrackAssoc const& trackIndices,
                      TracksPidWithSelAndMc const& tracks,
                      aod::McParticles const& particlesMc,
                      BCsInfo const& bcs,
                      McCollisions const& mcCollisions)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBplus(invMassWindowJpsiHad.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, DecayChannel::BplusToJpsiK>(collision, candsJpsiThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BplusToJpsiK>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorJpsiHadReduced, processJpsiKMc, "Process J/Psi K with MC info", false);

  void processJpsiPhiMc(CollisionsWCMcLabels const& collisions,
                        aod::HfCand2ProngWPid const& candsJpsi,
                        aod::TrackAssoc const& trackIndices,
                        TracksPidWithSelAndMc const& tracks,
                        aod::McParticles const& particlesMc,
                        BCsInfo const& bcs,
                        McCollisions const& mcCollisions)
  {
    // store configurables needed for B+ workflow
    if (!isHfCandBhadConfigFilled) {
      rowCandidateConfigBs(invMassWindowJpsiHad.value);
      isHfCandBhadConfigFilled = true;
    }

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);

      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto collsSameMcCollision = collisions.sliceBy(colPerMcCollision, collision.mcCollisionId());
      int64_t const indexCollisionMaxNumContrib = getIndexCollisionMaxNumContrib(collsSameMcCollision);
      runDataCreation<true, DecayChannel::BsToJpsiPhi>(collision, candsJpsiThisColl, trackIdsThisCollision, tracks, particlesMc, indexCollisionMaxNumContrib, bcs);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
    for (const auto& mcCollision : mcCollisions) {
      runMcGen<DecayChannel::BsToJpsiPhi>(mcCollision, particlesMc, collisions, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorJpsiHadReduced, processJpsiPhiMc, "Process J/Psi phi with MC info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorJpsiHadReduced>(cfgc)};
}
